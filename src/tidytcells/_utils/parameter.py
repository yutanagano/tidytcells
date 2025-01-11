import warnings


class Parameter:
    def __init__(self, value: any, name: str) -> None:
        self.value = value
        self.name = name

        if value is not None:
            self.is_explicitly_set = True
        else:
            self.is_explicitly_set = False

    def set_default(self, default_value: any) -> "Parameter":
        if not self.is_explicitly_set:
            self.value = default_value

        return self

    def resolve_with_alias(self, alias_value: any, alias_name: str) -> "Parameter":
        if alias_value is None:
            return self

        if self.is_explicitly_set:
            warnings.warn(
                f'The parameters "{self.name}" and "{alias_name}" are mutually exclusive. '
                "The latter will be deprecated in the near future, and replaced by the former. "
                "Both parameters were passed values, so disregarding that of the latter. "
                f'Please exclusively use "{self.name}" in the future.',
                FutureWarning,
            )
            return self

        warnings.warn(
            f'The parameter "{alias_name}" will be deprecated in the near future. Please use "{self.name}" instead.',
            FutureWarning,
        )
        self.value = alias_value

        return self

    def throw_error_if_not_of_type(
        self, object_type: type, optional: bool = False
    ) -> "Parameter":
        if optional and self.value is None:
            return self

        if not isinstance(self.value, object_type):
            raise TypeError(
                f'"{self.name}" must be of type {object_type}, got {self.value} ({type(self.value)}).'
            )

        return self

    def throw_error_if_not_one_of(self, *allowed_values) -> "Parameter":
        if not self.value in allowed_values:
            raise ValueError(
                f'"{self.name}" must be one of ({allowed_values}), got {self.value}.'
            )

        return self
