import re
import string
from typing import List

_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_PRODUCT_CODE_RE = re.compile(r"^(?=.{8,}$)(?=.*[A-Z])(?=.*\d)[A-Z0-9]+$")

FORBIDDEN_TERMS = ["caswell", "unii", "einecs", "ec:"]
FORBIDDEN_START_TERMS = ["ec ", "ec-", "cas ", "cas-", "wln ", "wln:"]


def cas_check_digit(potential_cas_number: str) -> bool:
    check_digit = potential_cas_number[-1]
    reversed_remaining_digits = potential_cas_number[::-1][1:]
    total = sum(
        (i + 1) * int(digit) for i, digit in enumerate(reversed_remaining_digits)
    )
    return total % 10 == int(check_digit)


def is_valid_cas(synonym: str) -> bool:
    if not isinstance(synonym, str):
        return False
    synonym = synonym.replace("-", "")
    if not synonym.isdigit():
        return False
    if len(synonym) < 3:
        return False
    return cas_check_digit(synonym)


def is_inchi_key(synonym: str, normalize: bool = True) -> bool:
    if not isinstance(synonym, str):
        return False
    if normalize:
        synonym = synonym.strip().upper()
    return _INCHIKEY_RE.match(synonym) is not None


def is_inchi(synonym: str) -> bool:
    return "inchi=" in synonym.lower()


def synonym_too_short(synonym: str, min_length: int = 3) -> bool:
    return len(synonym) < min_length


def synonym_too_long(synonym: str, max_length: int = 256) -> bool:
    return len(synonym) > max_length


def contains_percent(synonym: str) -> bool:
    return "%" in synonym


def has_only_numbers_punctuation_or_space(synonym: str) -> bool:
    return all(
        (ch.isdigit() or ch in string.punctuation or ch.isspace()) for ch in synonym
    )


def has_too_many_consecutive_digits(synonym: str, max_digits: int = 3) -> bool:
    digits = re.findall(r"\d+", synonym)
    if digits and len(max(digits, key=len)) > max_digits:
        return True
    return False


def looks_like_product_code(synonym: str, normalize: bool = True) -> bool:
    if normalize:
        synonym = synonym.strip().upper()
    return _PRODUCT_CODE_RE.fullmatch(synonym) is not None


def contains_forbidden_terms(
    synonym: str,
    forbidden_terms: List[str] = FORBIDDEN_TERMS,
    forbidden_start_terms: List[str] = FORBIDDEN_START_TERMS,
) -> bool:
    for forbidden_term in forbidden_terms:
        if forbidden_term.lower() in synonym.lower():
            return True
    for forbidden_start_term in forbidden_start_terms:
        if synonym.lower().startswith(forbidden_start_term.lower()):
            return True
    return False


def clean_synonyms_list(
    synonyms_list: List[str] | str,
    max_number_of_synonyms: int = 75,
) -> List[str]:
    if isinstance(synonyms_list, str):
        synonyms_list = [synonyms_list]
    if not synonyms_list:
        return []
    cleaned = []
    for synonym in synonyms_list:
        if synonym_too_short(synonym):
            continue
        if synonym_too_long(synonym):
            continue
        if is_valid_cas(synonym):
            continue
        if is_inchi_key(synonym):
            continue
        if is_inchi(synonym):
            continue
        if contains_percent(synonym):
            continue
        if has_only_numbers_punctuation_or_space(synonym):
            continue
        if has_too_many_consecutive_digits(synonym):
            continue
        if looks_like_product_code(synonym):
            continue
        if contains_forbidden_terms(synonym):
            continue
        cleaned.append(synonym)
    return list(set(cleaned))[:max_number_of_synonyms]


def get_cas_nos_from_synonyms_list(synonyms_list: List[str] | str) -> List[str]:
    if isinstance(synonyms_list, str):
        synonyms_list = [synonyms_list]
    if not synonyms_list:
        return []
    return list({s for s in synonyms_list if is_valid_cas(s)})
