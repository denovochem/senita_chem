"""Unit tests for senita_chem.synonym_cleaning module."""

from senita_chem.synonym_cleaning import (
    cas_check_digit,
    clean_synonyms_list,
    contains_forbidden_terms,
    contains_percent,
    get_cas_nos_from_synonyms_list,
    has_only_numbers_punctuation_or_space,
    has_too_many_consecutive_digits,
    is_inchi,
    is_inchi_key,
    is_valid_cas,
    looks_like_product_code,
    synonym_too_long,
    synonym_too_short,
)


class TestCasCheckDigit:
    """Tests for cas_check_digit function."""

    def test_valid_cas_ethanol(self) -> None:
        """cas_check_digit returns True for valid CAS 64-17-5."""
        assert cas_check_digit("64175") is True

    def test_valid_cas_longer(self) -> None:
        """cas_check_digit returns True for valid CAS with more digits."""
        assert cas_check_digit("7732185") is True

    def test_invalid_check_digit(self) -> None:
        """cas_check_digit returns False when check digit is wrong."""
        assert cas_check_digit("64170") is False

    def test_single_digit(self) -> None:
        """cas_check_digit returns True for a single digit (mod 10 matches)."""
        assert cas_check_digit("0") is True


class TestIsValidCas:
    """Tests for is_valid_cas function."""

    def test_valid_with_hyphens(self) -> None:
        """is_valid_cas validates CAS with hyphens."""
        assert is_valid_cas("64-17-5") is True

    def test_valid_without_hyphens(self) -> None:
        """is_valid_cas validates CAS without hyphens."""
        assert is_valid_cas("7732185") is True

    def test_invalid_check_digit(self) -> None:
        """is_valid_cas rejects CAS with bad check digit."""
        assert is_valid_cas("64-17-0") is False

    def test_non_numeric(self) -> None:
        """is_valid_cas rejects non-numeric input."""
        assert is_valid_cas("abc") is False

    def test_too_short(self) -> None:
        """is_valid_cas rejects strings shorter than 3 digits."""
        assert is_valid_cas("12") is False

    def test_non_string_input(self) -> None:
        """is_valid_cas returns False for non-string input."""
        assert is_valid_cas(12345) is False  # type: ignore[arg-type]

    def test_empty_string(self) -> None:
        """is_valid_cas rejects empty string."""
        assert is_valid_cas("") is False


class TestIsInchiKey:
    """Tests for is_inchi_key function."""

    def test_valid_inchikey(self) -> None:
        """is_inchi_key recognizes valid InChIKey."""
        assert is_inchi_key("LFQSCWFLJHTTHZ-UHFFFAOYSA-N") is True

    def test_valid_with_normalization(self) -> None:
        """is_inchi_key normalizes and recognizes InChIKey."""
        assert is_inchi_key("  lfqsCWFLJHTTHZ-uhfffaOYSA-n  ") is True

    def test_invalid_format(self) -> None:
        """is_inchi_key rejects malformed InChIKey."""
        assert is_inchi_key("NOT-AN-INCHIKEY") is False

    def test_short_string(self) -> None:
        """is_inchi_key rejects too-short string."""
        assert is_inchi_key("ABC") is False

    def test_non_string(self) -> None:
        """is_inchi_key returns False for non-string input."""
        assert is_inchi_key(12345) is False  # type: ignore[arg-type]

    def test_no_normalize(self) -> None:
        """is_inchi_key without normalization requires uppercase."""
        assert is_inchi_key("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", normalize=False) is True
        assert is_inchi_key("lfqscwfljhtthz-uhfffaoysa-n", normalize=False) is False


class TestIsInchi:
    """Tests for is_inchi function."""

    def test_valid_inchi(self) -> None:
        """is_inchi detects InChI prefix."""
        assert is_inchi("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3") is True

    def test_lowercase_prefix(self) -> None:
        """is_inchi is case-insensitive."""
        assert is_inchi("inchi=1/CH4/h1H4") is True

    def test_not_inchi(self) -> None:
        """is_inchi returns False for non-InChI string."""
        assert is_inchi("ethanol") is False

    def test_empty_string(self) -> None:
        """is_inchi returns False for empty string."""
        assert is_inchi("") is False


class TestSynonymTooShort:
    """Tests for synonym_too_short function."""

    def test_below_default(self) -> None:
        """synonym_too_short returns True for length < 3."""
        assert synonym_too_short("ab") is True

    def test_at_boundary(self) -> None:
        """synonym_too_short returns False for length == 3."""
        assert synonym_too_short("abc") is False

    def test_custom_min(self) -> None:
        """synonym_too_short respects custom min_length."""
        assert synonym_too_short("hello", min_length=10) is True


class TestSynonymTooLong:
    """Tests for synonym_too_long function."""

    def test_above_default(self) -> None:
        """synonym_too_long returns True for length > 256."""
        assert synonym_too_long("x" * 257) is True

    def test_at_boundary(self) -> None:
        """synonym_too_long returns False for length == 256."""
        assert synonym_too_long("x" * 256) is False

    def test_custom_max(self) -> None:
        """synonym_too_long respects custom max_length."""
        assert synonym_too_long("hello", max_length=4) is True


class TestContainsPercent:
    """Tests for contains_percent function."""

    def test_contains_percent(self) -> None:
        """contains_percent returns True when % is present."""
        assert contains_percent("50% solution") is True

    def test_no_percent(self) -> None:
        """contains_percent returns False when % is absent."""
        assert contains_percent("pure ethanol") is False


class TestHasOnlyNumbersPunctuationOrSpace:
    """Tests for has_only_numbers_punctuation_or_space function."""

    def test_only_digits(self) -> None:
        """has_only_numbers_punctuation_or_space returns True for only digits."""
        assert has_only_numbers_punctuation_or_space("1234") is True

    def test_digits_and_punctuation(self) -> None:
        """has_only_numbers_punctuation_or_space returns True for digits + punctuation."""
        assert has_only_numbers_punctuation_or_space("1,234.56") is True

    def test_with_letters(self) -> None:
        """has_only_numbers_punctuation_or_space returns False when letters present."""
        assert has_only_numbers_punctuation_or_space("1,234a") is False

    def test_empty_string(self) -> None:
        """has_only_numbers_punctuation_or_space returns True for empty string."""
        assert has_only_numbers_punctuation_or_space("") is True


class TestHasTooManyConsecutiveDigits:
    """Tests for has_too_many_consecutive_digits function."""

    def test_over_default(self) -> None:
        """has_too_many_consecutive_digits returns True for > 3 consecutive digits."""
        assert has_too_many_consecutive_digits("test1234") is True

    def test_at_default_boundary(self) -> None:
        """has_too_many_consecutive_digits returns False for exactly 3 digits."""
        assert has_too_many_consecutive_digits("test123") is False

    def test_custom_max(self) -> None:
        """has_too_many_consecutive_digits respects custom max_digits."""
        assert has_too_many_consecutive_digits("test12", max_digits=1) is True

    def test_no_digits(self) -> None:
        """has_too_many_consecutive_digits returns False when no digits."""
        assert has_too_many_consecutive_digits("nodigits") is False


class TestLooksLikeProductCode:
    """Tests for looks_like_product_code function."""

    def test_valid_product_code(self) -> None:
        """looks_like_product_code recognizes product codes."""
        assert looks_like_product_code("AB12345C") is True

    def test_too_short(self) -> None:
        """looks_like_product_code rejects short strings."""
        assert looks_like_product_code("AB12") is False

    def test_no_letters(self) -> None:
        """looks_like_product_code rejects all-digit strings."""
        assert looks_like_product_code("123456789") is False

    def test_no_digits(self) -> None:
        """looks_like_product_code rejects all-letter strings."""
        assert looks_like_product_code("ABCDEFGH") is False

    def test_lowercase_normalized(self) -> None:
        """looks_like_product_code normalizes to uppercase."""
        assert looks_like_product_code("ab12345c") is True

    def test_no_normalize(self) -> None:
        """looks_like_product_code without normalization is case-sensitive."""
        assert looks_like_product_code("ab12345c", normalize=False) is False


class TestContainsForbiddenTerms:
    """Tests for contains_forbidden_terms function."""

    def test_forbidden_mid(self) -> None:
        """contains_forbidden_terms detects mid-string forbidden terms."""
        assert contains_forbidden_terms("caswell no. 123") is True

    def test_forbidden_start(self) -> None:
        """contains_forbidden_terms detects start-of-string forbidden terms."""
        assert contains_forbidden_terms("ec number 123") is True

    def test_case_insensitive(self) -> None:
        """contains_forbidden_terms is case-insensitive."""
        assert contains_forbidden_terms("CASWELL NUMBER") is True

    def test_no_forbidden(self) -> None:
        """contains_forbidden_terms returns False for clean synonyms."""
        assert contains_forbidden_terms("ethyl alcohol") is False

    def test_custom_terms(self) -> None:
        """contains_forbidden_terms respects custom forbidden lists."""
        assert (
            contains_forbidden_terms(
                "badword here", forbidden_terms=["badword"], forbidden_start_terms=[]
            )
            is True
        )


class TestCleanSynonymsList:
    """Tests for clean_synonyms_list function."""

    def test_basic_filtering(self) -> None:
        """clean_synonyms_list removes CAS, InChI, short, and forbidden terms."""
        synonyms = [
            "ethanol",
            "64-17-5",
            "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
            "InChI=1S/C2H6O",
            "ab",
            "caswell no. 123",
            "ethyl alcohol",
            "50% solution",
            "AB12345C",
        ]
        result = clean_synonyms_list(synonyms)
        assert "ethanol" in result
        assert "ethyl alcohol" in result
        assert "64-17-5" not in result
        assert "LFQSCWFLJHTTHZ-UHFFFAOYSA-N" not in result
        assert "InChI=1S/C2H6O" not in result
        assert "ab" not in result
        assert "caswell no. 123" not in result
        assert "50% solution" not in result
        assert "AB12345C" not in result

    def test_string_input(self) -> None:
        """clean_synonyms_list handles single string input."""
        result = clean_synonyms_list("ethanol")
        assert result == ["ethanol"]

    def test_empty_list(self) -> None:
        """clean_synonyms_list returns empty list for empty input."""
        assert clean_synonyms_list([]) == []

    def test_deduplication(self) -> None:
        """clean_synonyms_list deduplicates synonyms."""
        result = clean_synonyms_list(["ethanol", "ethanol", "ethyl alcohol"])
        assert len(result) == 2
        assert "ethanol" in result
        assert "ethyl alcohol" in result

    def test_max_synonyms(self) -> None:
        """clean_synonyms_list respects max_number_of_synonyms."""
        synonyms = ["apple", "banana", "cherry", "date", "elderberry"]
        result = clean_synonyms_list(synonyms, max_number_of_synonyms=3)
        assert len(result) == 3

    def test_no_digits_only(self) -> None:
        """clean_synonyms_list removes digit-only/punctuation strings."""
        result = clean_synonyms_list(["ethanol", "1,2,3", "!!"])
        assert result == ["ethanol"]

    def test_too_long_rejected(self) -> None:
        """clean_synonyms_list removes overly long synonyms."""
        result = clean_synonyms_list(["ethanol", "x" * 257])
        assert result == ["ethanol"]

    def test_too_many_consecutive_digits(self) -> None:
        """clean_synonyms_list removes synonyms with >3 consecutive digits."""
        result = clean_synonyms_list(["ethanol", "version 1234"])
        assert result == ["ethanol"]


class TestGetCasNosFromSynonymsList:
    """Tests for get_cas_nos_from_synonyms_list function."""

    def test_extracts_cas(self) -> None:
        """get_cas_nos_from_synonyms_list extracts valid CAS numbers."""
        result = get_cas_nos_from_synonyms_list(["ethanol", "64-17-5", "ethyl alcohol"])
        assert result == ["64-17-5"]

    def test_string_input(self) -> None:
        """get_cas_nos_from_synonyms_list handles single string input."""
        result = get_cas_nos_from_synonyms_list("64-17-5")
        assert result == ["64-17-5"]

    def test_empty_list(self) -> None:
        """get_cas_nos_from_synonyms_list returns empty list for empty input."""
        assert get_cas_nos_from_synonyms_list([]) == []

    def test_deduplication(self) -> None:
        """get_cas_nos_from_synonyms_list deduplicates CAS numbers."""
        result = get_cas_nos_from_synonyms_list(["64-17-5", "64-17-5"])
        assert len(result) == 1

    def test_invalid_cas_excluded(self) -> None:
        """get_cas_nos_from_synonyms_list rejects invalid CAS numbers."""
        result = get_cas_nos_from_synonyms_list(["64-17-0", "ethanol"])
        assert result == []
