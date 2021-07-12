"""Customized PyGeoOGC exceptions."""
from typing import Generator, List, Optional, Union


class ServiceError(Exception):
    """Exception raised when the requested data is not available on the server.

    Parameters
    ----------
    err : str
        Service error message.
    """

    def __init__(self, err: str) -> None:
        self.message = f"Service returned the following error message:\n{err}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ServiceUnavailable(Exception):
    """Exception raised when the service is not available.

    Parameters
    ----------
    url : str
        The server url
    """

    def __init__(self, url: str) -> None:
        self.message = f"Service is currently not available, try again later:\n{url}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class ZeroMatched(ValueError):
    """Exception raised when a function argument is missing.

    Parameters
    ----------
    msg : str
        The exception error message
    """

    def __init__(self, msg: Optional[str] = None) -> None:
        if msg is None:
            self.message = "Service returned no features."
        else:
            self.message = f"Service returned no features with the following error message:\n{msg}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class InvalidInputValue(Exception):
    """Exception raised for invalid input.

    Parameters
    ----------
    inp : str
        Name of the input parameter
    valid_inputs : tuple
        List of valid inputs
    """

    def __init__(
        self, inp: str, valid_inputs: Union[List[str], Generator[str, None, None]]
    ) -> None:
        self.message = f"Given {inp} is invalid. Valid options are:\n" + "\n".join(
            str(i) for i in valid_inputs
        )
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class InvalidInputType(Exception):
    """Exception raised when a function argument type is invalid.

    Parameters
    ----------
    arg : str
        Name of the function argument
    valid_type : str
        The valid type of the argument
    example : str, optional
        An example of a valid form of the argument, defaults to None.
    """

    def __init__(self, arg: str, valid_type: str, example: Optional[str] = None) -> None:
        self.message = f"The {arg} argument should be of type {valid_type}"
        if example is not None:
            self.message += f":\n{example}"
        super().__init__(self.message)

    def __str__(self) -> str:
        return self.message


class MissingInputs(ValueError):
    """Exception raised when there are missing function arguments."""
