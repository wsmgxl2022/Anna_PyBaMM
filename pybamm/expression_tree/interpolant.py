#
# Interpolating class
#
import numpy as np
from scipy import interpolate
import warnings

import pybamm


class Interpolant(pybamm.Function):
    """
    Interpolate data in 1D.

    Parameters
    ----------
    x : iterable of :class:`numpy.ndarray`
        1-D array(s) of real values defining the data point coordinates.
    y : :class:`numpy.ndarray`
        The values of the function to interpolate at the data points.
    children : iterable of :class:`pybamm.Symbol`
        Node(s) to use when evaluating the interpolant. Each child corresponds to an
        entry of x
    name : str, optional
        Name of the interpolant. Default is None, in which case the name "interpolating
        function" is given.
    interpolator : str, optional
        Which interpolator to use. Can be "linear", "cubic", or "pchip". Default is
        "linear".
    extrapolate : bool, optional
        Whether to extrapolate for points that are outside of the parametrisation
        range, or return NaN (following default behaviour from scipy). Default is True.

    **Extends**: :class:`pybamm.Function`
    """

    def __init__(
        self,
        x,
        y,
        children,
        name=None,
        interpolator="linear",
        extrapolate=True,
        entries_string=None,
    ):
        # "cubic spline" has been renamed to "cubic"
        if interpolator == "cubic spline":
            interpolator = "cubic"
            warnings.warn(
                "The 'cubic spline' interpolator has been renamed to 'cubic'.",
                DeprecationWarning,
            )

        # Check interpolator is valid
        if interpolator not in ["linear", "cubic", "pchip"]:
            raise ValueError("interpolator '{}' not recognised".format(interpolator))

        # Perform some checks on the data
        if isinstance(x, (tuple, list)) and len(x) == 2:
            x1, x2 = x
            if y.ndim != 2:
                raise ValueError("y should be two-dimensional if len(x)=2")
            if x1.shape[0] != y.shape[1]:
                raise ValueError(
                    "len(x1) should equal y=shape[1], "
                    f"but x1.shape={x1.shape} and y.shape={y.shape}"
                )
            if x2 is not None and x2.shape[0] != y.shape[0]:
                raise ValueError(
                    "len(x2) should equal y=shape[0], "
                    f"but x2.shape={x2.shape} and y.shape={y.shape}"
                )
        else:
            if isinstance(x, (tuple, list)):
                x1 = x[0]
            else:
                x1 = x
                x = [x]
            x2 = None
            if x1.shape[0] != y.shape[0]:
                raise ValueError(
                    "len(x1) should equal y=shape[0], "
                    f"but x1.shape={x1.shape} and y.shape={y.shape}"
                )
        # children should be a list not a symbol
        if isinstance(children, pybamm.Symbol):
            children = [children]
        # Either a single x is provided and there is one child
        # or x is a 2-tuple and there are two children
        if len(x) != len(children):
            raise ValueError("len(x) should equal len(children)")
        # if there is only one x, y can be 2-dimensional but the child must have
        # length 1
        if len(x) == 1 and y.ndim == 2 and children[0].size != 1:
            raise ValueError(
                "child should have size 1 if y is two-dimensional and len(x)==1"
            )

        # Create interpolating function
        if len(x) == 1:
            self.dimension = 1
            if interpolator == "linear":
                if extrapolate is False:
                    fill_value = np.nan
                elif extrapolate is True:
                    fill_value = "extrapolate"
                interpolating_function = interpolate.interp1d(
                    x1,
                    y.T,
                    bounds_error=False,
                    fill_value=fill_value,
                )
            elif interpolator == "cubic":
                interpolating_function = interpolate.CubicSpline(
                    x1, y, extrapolate=extrapolate
                )
            elif interpolator == "pchip":
                interpolating_function = interpolate.PchipInterpolator(
                    x1, y, extrapolate=extrapolate
                )
        elif len(x) == 2:
            self.dimension = 2
            if interpolator == "pchip":
                raise ValueError(
                    "interpolator should be 'linear' or 'cubic' if x is two-dimensional"
                )
            else:
                interpolating_function = interpolate.interp2d(
                    x1, x2, y, kind=interpolator
                )
        else:
            raise ValueError("Invalid dimension of x: {0}".format(len(x)))

        # Set name
        if name is None:
            name = "interpolating_function"
        self.x = x
        self.y = y
        self.entries_string = entries_string
        super().__init__(
            interpolating_function, *children, name=name, derivative="derivative"
        )

        # Store information as attributes
        self.interpolator = interpolator
        self.extrapolate = extrapolate

    @property
    def entries_string(self):
        return self._entries_string

    @entries_string.setter
    def entries_string(self, value):
        # We must include the entries in the hash, since different arrays can be
        # indistinguishable by class, name and domain alone
        # Slightly different syntax for sparse and non-sparse matrices
        if value is not None:
            self._entries_string = value
        else:
            self._entries_string = ""
            for i, x in enumerate(self.x):
                self._entries_string += "x" + str(i) + "_" + str(x.tobytes())
            self._entries_string += "y_" + str(self.y.tobytes())

    def set_id(self):
        """See :meth:`pybamm.Symbol.set_id()`."""
        self._id = hash(
            (self.__class__, self.name, self.entries_string)
            + tuple([child.id for child in self.children])
            + tuple(self.domain)
        )

    def _function_new_copy(self, children):
        """See :meth:`Function._function_new_copy()`"""
        return pybamm.Interpolant(
            self.x,
            self.y,
            children,
            name=self.name,
            interpolator=self.interpolator,
            extrapolate=self.extrapolate,
            entries_string=self.entries_string,
        )

    def _function_evaluate(self, evaluated_children):
        children_eval_flat = []
        for child in evaluated_children:
            if isinstance(child, np.ndarray):
                children_eval_flat.append(child.flatten())
            else:
                children_eval_flat.append(child)
        if self.dimension == 1:
            return self.function(*children_eval_flat).flatten()[:, np.newaxis]
        elif self.dimension == 2:
            res = self.function(*children_eval_flat)
            if res.ndim > 1:
                return np.diagonal(res)[:, np.newaxis]
            else:
                # raise ValueError("Invalid children dimension: {0}".format(res.ndim))
                return res[:, np.newaxis]
        else:  # pragma: no cover
            raise ValueError("Invalid dimension: {0}".format(self.dimension))
