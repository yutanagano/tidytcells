import pytest
from tidytcells import aa


class TestStandardise:
    @pytest.mark.parametrize("seq", ("KLGGALQAK", "LLQTGIHVRVSQPSL", "SQLLNAKYL"))
    def test_already_correct(self, seq):
        result = aa.standardize(seq=seq)

        assert result == seq

    @pytest.mark.parametrize("seq", ("123456", "ASDFGHJKL", "A?AAAA", "AAAXAA"))
    def test_various_rejections(self, seq, caplog):
        result = aa.standardize(seq=seq)
        assert "not a valid amino acid" in caplog.text
        assert result == None

    @pytest.mark.parametrize(("seq", "expected"), (("klgak", "KLGAK"),))
    def test_various_corrections(self, seq, expected):
        result = aa.standardize(seq=seq)

        assert result == expected

    @pytest.mark.parametrize("seq", (None, 1, True, 5.4))
    def test_bad_input_type(self, seq):
        with pytest.raises(TypeError):
            aa.standardize(seq=seq)

    def test_standardise(self):
        result = aa.standardise(seq="KLGAK")

        assert result == "KLGAK"

    def test_log_failures(self, caplog):
        aa.standardize(seq="123456", log_failures=False)
        assert len(caplog.records) == 0

    def test_on_fail(self, caplog):
        result = aa.standardize("foobarbaz", on_fail="keep")
        assert "Failed to standardize" in caplog.text
        assert result == "foobarbaz"
