"""English prose policy: reject Han text without altering scientific identifiers."""


def contains_han(text):
    return any(0x3400 <= ord(c) <= 0x4DBF or 0x4E00 <= ord(c) <= 0x9FFF
               or 0xF900 <= ord(c) <= 0xFAFF or 0x20000 <= ord(c) <= 0x323AF
               for c in text)
