"""Validate work-block delivery and export complete member identity manifests."""
import argparse
import json
from pathlib import Path

from .boundary_pair_review import readonly
from .final_work_blocks import members, read, sealed, sha, validate


def property_members(blocks, properties, block_id):
    blocks, properties = Path(blocks).resolve(), Path(properties).resolve()
    parent = read(blocks / 'report.json'); child = read(properties / 'report.json')
    if Path(child['parent_blocks']).resolve() != blocks or child['sources'][str(blocks / 'report.json')] != sha(blocks / 'report.json'):
        raise ValueError('Property/parent artifacts disagree')
    if block_id == 'special-exhaustive':
        yield from members(blocks, block_id)
        return
    with readonly(properties / 'property_blocks.sqlite') as db:
        if not db.execute('SELECT 1 FROM block WHERE block_id=?', (block_id,)).fetchone():
            raise ValueError('Unknown property block')
        db.execute('ATTACH DATABASE ? AS source', ((Path(parent['model']) / 'blocks.sqlite').as_uri() + '?mode=ro',))
        db.execute('ATTACH DATABASE ? AS work', ((blocks / 'work_blocks.sqlite').as_uri() + '?mode=ro',))
        for cid, mid, original, target, state, rotation, leaf, parent_id in db.execute('''
            SELECT m.cid,p.mid,p.group_id,COALESCE(o.target_group,p.group_id),
            COALESCE(o.state,'unchanged_definite'),COALESCE(o.rotation,0),p.node,m.parent_block
            FROM membership m JOIN source.point p ON p.cid=m.cid LEFT JOIN work.override o ON o.cid=m.cid
            WHERE m.block_id=?''', (block_id,)):
            yield dict(cid=cid, mid=mid, original_group=original, logical_class=target, state=state,
                       rotation=rotation, old_leaf=leaf, block_id=block_id, parent_block=parent_id)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['validate', 'export'])
    parser.add_argument('--blocks', type=Path, required=True)
    parser.add_argument('--properties', type=Path)
    parser.add_argument('--block-id')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a fresh output file: ' + str(args.output))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.command == 'validate':
        result = validate(args.blocks, args.properties)
        args.output.write_text(json.dumps(result, indent=2), encoding='utf-8')
    else:
        if not args.block_id:
            raise ValueError('Export requires --block-id')
        br = sealed(args.blocks, ['work_blocks.sqlite'])
        sealed(br['model'], ['blocks.sqlite'])
        if args.properties:
            sealed(args.properties, ['property_blocks.sqlite'])
        iterable = property_members(args.blocks, args.properties, args.block_id) if args.properties else members(args.blocks, args.block_id)
        count = 0; temporary = args.output.with_name(args.output.name + '.partial')
        with temporary.open('x', encoding='utf-8') as stream:
            for row in iterable:
                stream.write(json.dumps(row) + '\n'); count += 1
        selected_database = args.properties / 'property_blocks.sqlite' if args.properties and args.block_id != 'special-exhaustive' else args.blocks / 'work_blocks.sqlite'
        with readonly(selected_database) as db:
            expected = db.execute('SELECT n FROM block WHERE block_id=?', (args.block_id,)).fetchone()[0]
        if count != expected:
            raise ValueError('Exported member count differs from declared block size')
        temporary.rename(args.output)
        result = dict(status='complete', conformers=count, output=str(args.output), sha256=sha(args.output))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
