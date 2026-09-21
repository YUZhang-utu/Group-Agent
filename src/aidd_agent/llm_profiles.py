"""Resolve one explicitly selected provider; never silently fail over."""
import os
from pathlib import Path


def select_llm_profile(provider=None, explicit=None, config_dir=None):
    if provider is not None and provider not in {"gpt", "deepseek"}:
        raise ValueError("Unknown provider")
    if provider and explicit:
        raise ValueError("Choose provider OR explicit LLM profile")
    if explicit:
        result = Path(explicit)
    elif provider is None and os.environ.get("AIDD_LLM_PROFILE"):
        result = Path(os.environ["AIDD_LLM_PROFILE"])
    else:
        name = provider or os.environ.get("AIDD_LLM_PROVIDER", "gpt")
        if name not in {"gpt", "deepseek"}:
            raise ValueError("Unknown provider")
        directory = config_dir or os.environ.get("AIDD_LLM_CONFIG_DIR")
        result = (Path(directory) if directory else Path(__file__).resolve().parents[2] / "configs/llm") / (name + ".json")
    if not result.is_file():
        raise ValueError(f"LLM profile not found: {result}")
    return result.resolve()
