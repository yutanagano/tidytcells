import warnings


class Parameter:
    def __init__(self, value: any, name: str) -> None:
        self.value = value
        self.name = name

    @property
    def is_specified(self) -> bool:
        return self.value is not None

    def resolve_with_alias_and_return_value(self, alias: "Parameter") -> any:
        if self.is_specified:
            return self.value

        if alias.is_specified:
            warnings.warn(
                f'The parameter "{alias.name}" will be deprecated in the near future. Please use "{self.name}" instead.',
                FutureWarning,
            )
            return alias.value

        return None

    def throw_error_if_not_of_type(
        self, object_type: type, optional: bool = False
    ) -> None:
        if optional and self.value is None:
            return

        if not isinstance(self.value, object_type):
            raise TypeError(
                f'"{self.name}" must be of type {object_type}, got {self.value} ({type(self.value)}).'
            )

    def throw_error_if_not_one_of(self, *allowed_values) -> None:
        if not self.value in allowed_values:
            raise ValueError(
                f'"{self.name}" must be one of ({allowed_values}), got {self.value}.'
            )
