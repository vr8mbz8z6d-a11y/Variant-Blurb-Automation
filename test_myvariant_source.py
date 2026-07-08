import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from myvariant_source import parse_gnomad, parse_revel, fetch_myvariant_record
from unittest.mock import patch, MagicMock

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "..", "fixtures", "myvariant_sample.json")


def test_fetch_myvariant_record_includes_assembly_param_for_hg38():
    """
    Locks in the real bug from the live run: MyVariant.info defaults to
    hg19 unless assembly=hg38 is explicitly passed. Without this, an hg38
    coordinate gets silently looked up against hg19's map and returns
    nothing -- exactly what happened the first time this pipeline ran
    against a real variant.
    """
    fake_resp = MagicMock(status_code=200, text="{}")
    fake_resp.json.return_value = {}
    with patch("myvariant_source.requests.get", return_value=fake_resp) as mock_get:
        fetch_myvariant_record("chr13:g.32370447G>A", genome_build="hg38")
        _, kwargs = mock_get.call_args
        assert kwargs["params"].get("assembly") == "hg38"


def test_fetch_myvariant_record_omits_assembly_param_for_hg19():
    fake_resp = MagicMock(status_code=200, text="{}")
    fake_resp.json.return_value = {}
    with patch("myvariant_source.requests.get", return_value=fake_resp) as mock_get:
        fetch_myvariant_record("chr13:g.32944584G>A", genome_build="hg19")
        _, kwargs = mock_get.call_args
        assert "assembly" not in kwargs["params"]


def load_fixture():
    with open(FIXTURE_PATH) as f:
        return json.load(f)


def test_parse_gnomad():
    record = load_fixture()
    g = parse_gnomad(record)
    assert g.af == 0.0000124
    assert g.allele_count == 3
    assert g.allele_number == 241830
    assert g.homozygote_count == 0


def test_parse_revel():
    record = load_fixture()
    r = parse_revel(record)
    assert r.score == 0.9


def test_parse_gnomad_missing_data_returns_none_fields():
    g = parse_gnomad({})
    assert g.af is None
    assert g.allele_count is None


def test_parse_revel_missing_data_returns_none():
    r = parse_revel({})
    assert r.score is None


def test_parse_revel_handles_list_of_scores():
    r = parse_revel({"dbnsfp": {"revel": {"score": [0.4, 0.9, 0.6]}}})
    assert r.score == 0.9  # takes the max, the conservative read
