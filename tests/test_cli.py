"""Tests for the summarize-index CLI entry point."""

import json
from pathlib import Path

import pytest

from ngsindex import cli

DATA_DIR = Path(__file__).parent / "test_index"


def test_summarize_to_file(tmp_path, monkeypatch):
    out = tmp_path / "summary.json"
    monkeypatch.setattr(
        "sys.argv",
        ["summarize-index", "-o", str(out), str(DATA_DIR / "small.bam.bai")],
    )
    cli.summarize()
    summary = json.loads(out.read_text())
    assert "num_references" in summary
    assert "read_counts" in summary


def test_summarize_to_stdout(capsys, monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        ["summarize-index", str(DATA_DIR / "small.vcf.gz.tbi")],
    )
    cli.summarize()
    captured = capsys.readouterr()
    summary = json.loads(captured.out)
    assert summary["num_references"] >= 1


def test_summarize_requires_index_arg(monkeypatch):
    monkeypatch.setattr("sys.argv", ["summarize-index"])
    with pytest.raises(SystemExit):
        cli.summarize()
