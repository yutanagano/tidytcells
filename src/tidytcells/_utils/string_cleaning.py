def clean_and_lowercase(string: str) -> str:
    string = remove_whitespace_and_pollutants(string)
    return string.lower()


def clean_and_uppercase(string: str) -> str:
    string = remove_whitespace_and_pollutants(string)
    return string.upper()


def remove_whitespace_and_pollutants(string: str) -> str:
    string = "".join(string.split())
    string = string.replace("&nbsp;", "")
    string = string.replace("&ndash;", "-")
    return string
