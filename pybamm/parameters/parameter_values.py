#
# Dimensional and dimensionless parameter values, and scales
#
import numpy as np
import pybamm
import pandas as pd
import os
import numbers
import warnings
from pprint import pformat
from collections import defaultdict
import inspect
import json


class ParameterValues:
    """
    The parameter values for a simulation.

    Note that this class does not inherit directly from the python dictionary class as
    this causes issues with saving and loading simulations.

    Parameters
    ----------
    values : dict or string
        Explicit set of parameters, or reference to a file of parameters
        If string, gets passed to read_parameters_csv to read a file.
    chemistry : dict
        Dict of strings for default chemistries. Must be of the form:
        {"base chemistry": base_chemistry,
        "cell": cell_properties_authorYear,
        "negative electrode": negative_electrode_chemistry_authorYear,
        "separator": separator_chemistry_authorYear,
        "positive electrode": positive_electrode_chemistry_authorYear,
        "electrolyte": electrolyte_chemistry_authorYear,
        "experiment": experimental_conditions_authorYear}.
        Then the negative electrode chemistry is loaded from the file
        inputs/parameters/base_chemistry/negative electrodes/
        negative_electrode_chemistry_authorYear, etc.
        Parameters in "cell" should include geometry and current collector properties.
        Parameters in "experiment" should include parameters relating to experimental
        conditions, such as initial conditions and currents.

    Examples
    --------
    >>> import pybamm
    >>> values = {"some parameter": 1, "another parameter": 2}
    >>> param = pybamm.ParameterValues(values)
    >>> param["some parameter"]
    1
    >>> file = "input/parameters/lithium_ion/cells/kokam_Marquis2019/parameters.csv"
    >>> values_path = pybamm.get_parameters_filepath(file)
    >>> param = pybamm.ParameterValues(values_path)
    >>> param["Negative current collector thickness [m]"]
    2.5e-05
    >>> param = pybamm.ParameterValues("Marquis2019")
    >>> param["Reference temperature [K]"]
    298.15

    """

    def __init__(self, values=None, chemistry=None):
        self._dict_items = pybamm.FuzzyDict()
        # Must provide either values or chemistry, not both (nor neither)
        if values is not None and chemistry is not None:
            raise ValueError(
                "Only one of values and chemistry can be provided. To change parameters"
                " slightly from a chemistry, first load parameters with the chemistry"
                " (param = pybamm.ParameterValues(...)) and then update with"
                " param.update({dict of values})."
            )
        if values is None and chemistry is None:
            raise ValueError("values and chemistry cannot both be None")
        # First load chemistry
        if chemistry is not None:
            warnings.warn(
                "The 'chemistry' keyword argument has been deprecated and will be "
                "removed in a future release. Call `ParameterValues` with a "
                "parameter set dictionary, or the name of a parameter set (string), "
                "as the single argument, e.g. `ParameterValues('Chen2020')`.",
                DeprecationWarning,
            )
            self.update_from_chemistry(chemistry)
        # Then update with values dictionary or file
        if values is not None:
            if (isinstance(values, str) and hasattr(pybamm.parameter_sets, values)) or (
                isinstance(values, dict) and "chemistry" in values
            ):
                self.update_from_chemistry(values)
            else:
                # If base_parameters is a filename, load from that filename
                if isinstance(values, str):
                    file_path = self.find_parameter(values)
                    path = os.path.split(file_path)[0]
                    values = self.read_parameters_csv(file_path)
                else:
                    path = ""
                # Don't check parameter already exists when first creating it
                self.update(values, check_already_exists=False, path=path)

        # Initialise empty _processed_symbols dict (for caching)
        self._processed_symbols = {}
        self.parameter_events = []

        # Don't touch this parameter unless you know what you are doing
        # This is for the conversion to Julia (ModelingToolkit)
        self._replace_callable_function_parameters = True

    def __getitem__(self, key):
        return self._dict_items[key]

    def get(self, key, default=None):
        """Return item correspoonding to key if it exists, otherwise return default"""
        try:
            return self._dict_items[key]
        except KeyError:
            return default

    def __setitem__(self, key, value):
        """Call the update functionality when doing a setitem"""
        self.update({key: value})

    def __delitem__(self, key):
        del self._dict_items[key]

    def __repr__(self):
        return pformat(self._dict_items, width=1)

    def __eq__(self, other):
        return self._dict_items == other._dict_items

    def keys(self):
        """Get the keys of the dictionary"""
        return self._dict_items.keys()

    def values(self):
        """Get the values of the dictionary"""
        return self._dict_items.values()

    def items(self):
        """Get the items of the dictionary"""
        return self._dict_items.items()

    def copy(self):
        """Returns a copy of the parameter values. Makes sure to copy the internal
        dictionary."""
        new_copy = ParameterValues(self._dict_items.copy())
        new_copy._replace_callable_function_parameters = (
            self._replace_callable_function_parameters
        )
        return new_copy

    def search(self, key, print_values=True):
        """
        Search dictionary for keys containing 'key'.

        See :meth:`pybamm.FuzzyDict.search()`.
        """
        return self._dict_items.search(key, print_values)

    def update_from_chemistry(self, chemistry):
        """
        Load standard set of components from a 'chemistry' dictionary
        """
        if isinstance(chemistry, str):
            chemistry = getattr(pybamm.parameter_sets, chemistry)

        base_chemistry = chemistry["chemistry"]

        # Load each component name

        component_groups = [
            "cell",
            "negative electrode",
            "positive electrode",
            "separator",
            "electrolyte",
            "experiment",
        ]

        # add SEI parameters if provided
        for extra_group in ["sei", "lithium plating"]:
            if extra_group in chemistry:
                component_groups = [extra_group] + component_groups

        for component_group in component_groups:
            # Make sure component is provided
            try:
                component = chemistry[component_group]
            except KeyError:
                raise KeyError(
                    "must provide '{}' parameters for {} chemistry".format(
                        component_group, base_chemistry
                    )
                )
            # Create path to component and load values
            component_path = os.path.join(
                base_chemistry, component_group.replace(" ", "_") + "s", component
            )
            file_path = self.find_parameter(
                os.path.join(component_path, "parameters.csv")
            )
            component_params_tmp = self.read_parameters_csv(file_path)

            component_params = {}
            for k, v in component_params_tmp.items():
                component_params[k] = v

            # Update parameters, making sure to check any conflicts
            self.update(
                component_params,
                check_conflict=True,
                check_already_exists=False,
                path=os.path.dirname(file_path),
            )

        # register (list of) citations
        if "citation" in chemistry:
            citations = chemistry["citation"]
            if not isinstance(citations, list):
                citations = [citations]
            for citation in citations:
                pybamm.citations.register(citation)

    def read_parameters_csv(self, filename):
        """Reads parameters from csv file into dict.

        Parameters
        ----------
        filename : str
            The name of the csv file containing the parameters.

        Returns
        -------
        dict
            {name: value} pairs for the parameters.

        """
        df = pd.read_csv(filename, comment="#", skip_blank_lines=True)
        # Drop rows that are all NaN (seems to not work with skip_blank_lines)
        df.dropna(how="all", inplace=True)
        return {k: v for (k, v) in zip(df["Name [units]"], df["Value"])}

    def update(self, values, check_conflict=False, check_already_exists=True, path=""):
        """
        Update parameter dictionary, while also performing some basic checks.

        Parameters
        ----------
        values : dict
            Dictionary of parameter values to update parameter dictionary with
        check_conflict : bool, optional
            Whether to check that a parameter in `values` has not already been defined
            in the parameter class when updating it, and if so that its value does not
            change. This is set to True during initialisation, when parameters are
            combined from different sources, and is False by default otherwise
        check_already_exists : bool, optional
            Whether to check that a parameter in `values` already exists when trying to
            update it. This is to avoid cases where an intended change in the parameters
            is ignored due a typo in the parameter name, and is True by default but can
            be manually overridden.
        path : string, optional
            Path from which to load functions
        """
        # check if values is not a dictionary
        if not isinstance(values, dict):
            values = values._dict_items
        # check parameter values
        self.check_parameter_values(values)
        # update
        for name, value in values.items():
            # check for conflicts
            if (
                check_conflict is True
                and name in self.keys()
                and not (self[name] == float(value) or self[name] == value)
            ):
                raise ValueError(
                    "parameter '{}' already defined with value '{}'".format(
                        name, self[name]
                    )
                )
            # check parameter already exists (for updating parameters)
            if check_already_exists is True:
                try:
                    self._dict_items[name]
                except KeyError as err:
                    raise KeyError(
                        "Cannot update parameter '{}' as it does not ".format(name)
                        + "have a default value. ({}). If you are ".format(err.args[0])
                        + "sure you want to update this parameter, use "
                        + "param.update({{name: value}}, check_already_exists=False)"
                    )
            # if no conflicts, update, loading functions and data if they are specified
            # Functions are flagged with the string "[function]"
            if isinstance(value, str):
                if value.startswith("[function]"):
                    loaded_value = pybamm.load_function(os.path.join(path, value[10:]))
                    self._dict_items[name] = loaded_value
                # Data is flagged with the string "[data]" or "[current data]"
                elif value.startswith("[current data]") or value.startswith("[data]"):
                    if value.startswith("[current data]"):
                        data_path = os.path.join(
                            pybamm.root_dir(), "pybamm", "input", "drive_cycles"
                        )
                        filename = os.path.join(data_path, value[14:] + ".csv")
                        function_name = value[14:]
                    else:
                        filename = os.path.join(path, value[6:] + ".csv")
                        function_name = value[6:]
                    filename = pybamm.get_parameters_filepath(filename)
                    data = pd.read_csv(
                        filename, comment="#", skip_blank_lines=True, header=None
                    ).to_numpy()
                    # Save name and data
                    self._dict_items[name] = (function_name, ([data[:, 0]], data[:, 1]))
                # parse 2D parameter data
                elif value.startswith("[2D data]"):
                    filename = os.path.join(path, value[9:] + ".json")
                    function_name = value[9:]
                    filename = pybamm.get_parameters_filepath(filename)
                    with open(filename, "r") as jsonfile:
                        json_data = json.load(jsonfile)
                    data = json_data["data"]
                    data[0] = [np.array(el) for el in data[0]]
                    data[1] = np.array(data[1])
                    self._dict_items[name] = (function_name, data)

                elif value == "[input]":
                    self._dict_items[name] = pybamm.InputParameter(name)
                # Anything else should be a converted to a float
                else:
                    self._dict_items[name] = float(value)
            elif isinstance(value, tuple) and isinstance(value[1], np.ndarray):
                # If data is provided as a 2-column array (1D data),
                # convert to two arrays for compatibility with 2D data
                # see #1805
                func_name, data = value
                data = ([data[:, 0]], data[:, 1])
                self._dict_items[name] = (func_name, data)
            else:
                self._dict_items[name] = value
        # reset processed symbols
        self._processed_symbols = {}

    def check_parameter_values(self, values):
        # Make sure typical current is non-zero
        if "Typical current [A]" in values and values["Typical current [A]"] == 0:
            raise ValueError(
                "'Typical current [A]' cannot be zero. A possible alternative is to "
                "set 'Current function [A]' to `0` instead."
            )

        for param in values:
            if "propotional term" in param:
                raise ValueError(
                    f"The parameter '{param}' has been renamed to "
                    "'... proportional term [s-1]', and its value should now be divided"
                    "by 3600 to get the same results as before."
                )

    def process_model(self, unprocessed_model, inplace=True):
        """Assign parameter values to a model.
        Currently inplace, could be changed to return a new model.

        Parameters
        ----------
        unprocessed_model : :class:`pybamm.BaseModel`
            Model to assign parameter values for
        inplace: bool, optional
            If True, replace the parameters in the model in place. Otherwise, return a
            new model with parameter values set. Default is True.

        Raises
        ------
        :class:`pybamm.ModelError`
            If an empty model is passed (`model.rhs = {}` and `model.algebraic = {}` and
            `model.variables = {}`)

        """
        pybamm.logger.info(
            "Start setting parameters for {}".format(unprocessed_model.name)
        )

        # set up inplace vs not inplace
        if inplace:
            # any changes to unprocessed_model attributes will change model attributes
            # since they point to the same object
            model = unprocessed_model
        else:
            # create a copy of the model
            model = unprocessed_model.new_copy()

        if (
            len(unprocessed_model.rhs) == 0
            and len(unprocessed_model.algebraic) == 0
            and len(unprocessed_model.variables) == 0
        ):
            raise pybamm.ModelError("Cannot process parameters for empty model")

        new_rhs = {}
        for variable, equation in unprocessed_model.rhs.items():
            pybamm.logger.verbose(
                "Processing parameters for {!r} (rhs)".format(variable)
            )
            new_rhs[variable] = self.process_symbol(equation)
        model.rhs = new_rhs

        new_algebraic = {}
        for variable, equation in unprocessed_model.algebraic.items():
            pybamm.logger.verbose(
                "Processing parameters for {!r} (algebraic)".format(variable)
            )
            new_algebraic[variable] = self.process_symbol(equation)
        model.algebraic = new_algebraic

        new_initial_conditions = {}
        for variable, equation in unprocessed_model.initial_conditions.items():
            pybamm.logger.verbose(
                "Processing parameters for {!r} (initial conditions)".format(variable)
            )
            new_initial_conditions[variable] = self.process_symbol(equation)
        model.initial_conditions = new_initial_conditions

        model.boundary_conditions = self.process_boundary_conditions(unprocessed_model)

        new_variables = {}
        for variable, equation in unprocessed_model.variables.items():
            pybamm.logger.verbose(
                "Processing parameters for {!r} (variables)".format(variable)
            )
            new_variables[variable] = self.process_symbol(equation)
        model.variables = new_variables

        new_events = []
        for event in unprocessed_model.events:
            pybamm.logger.verbose(
                "Processing parameters for event '{}''".format(event.name)
            )
            new_events.append(
                pybamm.Event(
                    event.name, self.process_symbol(event.expression), event.event_type
                )
            )

        for event in self.parameter_events:
            pybamm.logger.verbose(
                "Processing parameters for event '{}''".format(event.name)
            )
            new_events.append(
                pybamm.Event(
                    event.name, self.process_symbol(event.expression), event.event_type
                )
            )

        model.events = new_events

        # Set external variables
        model.external_variables = [
            self.process_symbol(var) for var in unprocessed_model.external_variables
        ]

        # Process timescale
        new_timescale = self.process_symbol(unprocessed_model.timescale)
        if isinstance(new_timescale, pybamm.Scalar):
            model._timescale = new_timescale
        else:
            raise ValueError(
                "model.timescale must be a Scalar after parameter processing "
                "(cannot contain 'InputParameter's). "
                "You have probably set one of the parameters used to calculate the "
                "timescale to an InputParameter. To avoid this error, hardcode "
                "model.timescale to a constant value by passing the option "
                "{'timescale': value} to the model."
            )

        # Process length scales
        new_length_scales = {}
        for domain, scale in unprocessed_model.length_scales.items():
            new_scale = self.process_symbol(scale)
            new_length_scales[domain] = new_scale
        model._length_scales = new_length_scales

        pybamm.logger.info("Finish setting parameters for {}".format(model.name))

        return model

    def process_boundary_conditions(self, model):
        """
        Process boundary conditions for a model
        Boundary conditions are dictionaries {"left": left bc, "right": right bc}
        in general, but may be imposed on the tabs (or *not* on the tab) for a
        small number of variables, e.g. {"negative tab": neg. tab bc,
        "positive tab": pos. tab bc "no tab": no tab bc}.
        """
        new_boundary_conditions = {}
        sides = ["left", "right", "negative tab", "positive tab", "no tab"]
        for variable, bcs in model.boundary_conditions.items():
            processed_variable = self.process_symbol(variable)
            new_boundary_conditions[processed_variable] = {}
            for side in sides:
                try:
                    bc, typ = bcs[side]
                    pybamm.logger.verbose(
                        "Processing parameters for {!r} ({} bc)".format(variable, side)
                    )
                    processed_bc = (self.process_symbol(bc), typ)
                    new_boundary_conditions[processed_variable][side] = processed_bc
                except KeyError as err:
                    # don't raise error if the key error comes from the side not being
                    # found
                    if err.args[0] in side:
                        pass
                    # do raise error otherwise (e.g. can't process symbol)
                    else:
                        raise KeyError(err)

        return new_boundary_conditions

    def process_geometry(self, geometry):
        """
        Assign parameter values to a geometry (inplace).

        Parameters
        ----------
        geometry : dict
            Geometry specs to assign parameter values to
        """

        def process_and_check(sym):
            if isinstance(sym, numbers.Number):
                return pybamm.Scalar(sym)
            new_sym = self.process_symbol(sym)
            if not isinstance(new_sym, pybamm.Scalar):
                raise ValueError(
                    "Geometry parameters must be Scalars after parameter processing"
                )
            return new_sym

        for domain in geometry:
            for spatial_variable, spatial_limits in geometry[domain].items():
                # process tab information if using 1 or 2D current collectors
                if spatial_variable == "tabs":
                    for tab, position_size in spatial_limits.items():
                        for position_size, sym in position_size.items():
                            geometry[domain]["tabs"][tab][
                                position_size
                            ] = process_and_check(sym)
                else:
                    for lim, sym in spatial_limits.items():
                        geometry[domain][spatial_variable][lim] = process_and_check(sym)

    def process_symbol(self, symbol):
        """Walk through the symbol and replace any Parameter with a Value.
        If a symbol has already been processed, the stored value is returned.

        Parameters
        ----------
        symbol : :class:`pybamm.Symbol`
            Symbol or Expression tree to set parameters for

        Returns
        -------
        symbol : :class:`pybamm.Symbol`
            Symbol with Parameter instances replaced by Value

        """
        try:
            return self._processed_symbols[symbol]
        except KeyError:
            processed_symbol = self._process_symbol(symbol)
            self._processed_symbols[symbol] = processed_symbol

            return processed_symbol

    def _process_symbol(self, symbol):
        """See :meth:`ParameterValues.process_symbol()`."""

        if isinstance(symbol, pybamm.Parameter):
            value = self[symbol.name]
            if isinstance(value, numbers.Number):
                # Check not NaN (parameter in csv file but no value given)
                if np.isnan(value):
                    raise ValueError(f"Parameter '{symbol.name}' not found")
                # Scalar inherits name (for updating parameters)
                return pybamm.Scalar(value, name=symbol.name)
            elif isinstance(value, pybamm.Symbol):
                new_value = self.process_symbol(value)
                new_value.copy_domains(symbol)
                return new_value
            else:
                raise TypeError("Cannot process parameter '{}'".format(value))

        elif isinstance(symbol, pybamm.FunctionParameter):
            new_children = []
            for child in symbol.children:
                if symbol.diff_variable is not None and any(
                    x == symbol.diff_variable for x in child.pre_order()
                ):
                    # Wrap with NotConstant to avoid simplification,
                    # which would stop symbolic diff from working properly
                    new_child = pybamm.NotConstant(child)
                    new_children.append(self.process_symbol(new_child))
                else:
                    new_children.append(self.process_symbol(child))
            function_name = self[symbol.name]

            # Create Function or Interpolant or Scalar object
            if isinstance(function_name, tuple):
                if len(function_name) == 2:  # CSV or JSON parsed data
                    # to create an Interpolant
                    name, data = function_name

                    if len(data[0]) == 1:
                        input_data = data[0][0], data[1]

                    else:
                        input_data = data

                    # For parameters provided as data we use a cubic interpolant
                    # Note: the cubic interpolant can be differentiated
                    function = pybamm.Interpolant(
                        input_data[0],
                        input_data[-1],
                        new_children,
                        interpolator="cubic",
                        name=name,
                    )
                    # Define event to catch extrapolation. In these events the sign is
                    # important: it should be positive inside of the range and negative
                    # outside of it
                    for data_index in range(len(data[0])):
                        self.parameter_events.append(
                            pybamm.Event(
                                "Interpolant {} lower bound".format(name),
                                pybamm.min(
                                    new_children[data_index] - min(data[0][data_index])
                                ),
                                pybamm.EventType.INTERPOLANT_EXTRAPOLATION,
                            )
                        )
                        self.parameter_events.append(
                            pybamm.Event(
                                "Interpolant {} upper bound".format(name),
                                pybamm.min(
                                    max(data[0][data_index]) - new_children[data_index]
                                ),
                                pybamm.EventType.INTERPOLANT_EXTRAPOLATION,
                            )
                        )

                else:  # pragma: no cover
                    raise ValueError(
                        "Invalid function name length: {0}".format(len(function_name))
                    )

            elif isinstance(function_name, numbers.Number):
                # Check not NaN (parameter in csv file but no value given)
                if np.isnan(function_name):
                    raise ValueError(
                        f"Parameter '{symbol.name}' (possibly a function) not found"
                    )
                # If the "function" is provided is actually a scalar, return a Scalar
                # object instead of throwing an error.
                function = pybamm.Scalar(function_name, name=symbol.name)
            elif callable(function_name):
                # otherwise evaluate the function to create a new PyBaMM object
                function = function_name(*new_children)
                if (
                    self._replace_callable_function_parameters is False
                    and not isinstance(
                        self.process_symbol(function), (pybamm.Scalar, pybamm.Broadcast)
                    )
                    and symbol.print_name is not None
                    and symbol.diff_variable is None
                ):
                    # Special trick for printing in Julia ModelingToolkit format
                    out = pybamm.FunctionParameter(
                        symbol.print_name, dict(zip(symbol.input_names, new_children))
                    )

                    out.arg_names = inspect.getfullargspec(function_name)[0]
                    out.callable = self.process_symbol(
                        function_name(
                            *[
                                pybamm.Variable(arg_name, domains=child.domains)
                                for arg_name, child in zip(out.arg_names, new_children)
                            ]
                        )
                    )

                    return out
            elif isinstance(
                function_name, (pybamm.Interpolant, pybamm.InputParameter)
            ) or (
                isinstance(function_name, pybamm.Symbol)
                and function_name.size_for_testing == 1
            ):
                function = function_name
            else:
                raise TypeError(
                    "Parameter provided for '{}' ".format(symbol.name)
                    + "is of the wrong type (should either be scalar-like or callable)"
                )
            # Differentiate if necessary
            if symbol.diff_variable is None:
                # Use ones_like so that we get the right shapes
                function_out = function * pybamm.ones_like(*new_children)
            else:
                # return differentiated function
                new_diff_variable = self.process_symbol(symbol.diff_variable)
                function_out = function.diff(new_diff_variable)
            # Process again just to be sure
            return self.process_symbol(function_out)

        elif isinstance(symbol, pybamm.BinaryOperator):
            # process children
            new_left = self.process_symbol(symbol.left)
            new_right = self.process_symbol(symbol.right)
            # make new symbol, ensure domain remains the same
            new_symbol = symbol._binary_new_copy(new_left, new_right)
            new_symbol.copy_domains(symbol)
            return new_symbol

        # Unary operators
        elif isinstance(symbol, pybamm.UnaryOperator):
            new_child = self.process_symbol(symbol.child)
            new_symbol = symbol._unary_new_copy(new_child)
            # ensure domain remains the same
            new_symbol.copy_domains(symbol)
            # x_average can sometimes create a new symbol with electrode thickness
            # parameters, so we process again to make sure these parameters are set
            if isinstance(symbol, pybamm.XAverage) and not isinstance(
                new_symbol, pybamm.XAverage
            ):
                new_symbol = self.process_symbol(new_symbol)
            # f_a_dist in the size average needs to be processed
            if isinstance(new_symbol, pybamm.SizeAverage):
                new_symbol.f_a_dist = self.process_symbol(new_symbol.f_a_dist)
            return new_symbol

        # Functions
        elif isinstance(symbol, pybamm.Function):
            new_children = [self.process_symbol(child) for child in symbol.children]
            return symbol._function_new_copy(new_children)

        # Concatenations
        elif isinstance(symbol, pybamm.Concatenation):
            new_children = [self.process_symbol(child) for child in symbol.children]
            return symbol._concatenation_new_copy(new_children)

        else:
            # Backup option: return the object
            return symbol

    def evaluate(self, symbol):
        """
        Process and evaluate a symbol.

        Parameters
        ----------
        symbol : :class:`pybamm.Symbol`
            Symbol or Expression tree to evaluate

        Returns
        -------
        number or array
            The evaluated symbol
        """
        processed_symbol = self.process_symbol(symbol)
        if processed_symbol.is_constant():
            return processed_symbol.evaluate()
        else:
            raise ValueError("symbol must evaluate to a constant scalar or array")

    def _ipython_key_completions_(self):
        return list(self._dict_items.keys())

    def export_csv(self, filename):

        # process functions and data to output
        # like they appear in inputs csv files
        parameter_output = {}
        for key, val in self.items():
            if callable(val):
                val = "[function]" + val.__name__
            elif isinstance(val, tuple):
                val = "[data]" + val[0]
            parameter_output[key] = [val]

        df = pd.DataFrame(parameter_output)
        df = df.transpose()
        df.to_csv(filename, header=["Value"], index_label="Name [units]")

    def print_parameters(self, parameters, output_file=None):
        """
        Return dictionary of evaluated parameters, and optionally print these evaluated
        parameters to an output file.
        For dimensionless parameters that depend on the C-rate, the value is given as a
        function of the C-rate (either x * Crate or x / Crate depending on the
        dependence)

        Parameters
        ----------
        parameters : class or dict containing :class:`pybamm.Parameter` objects
            Class or dictionary containing all the parameters to be evaluated
        output_file : string, optional
            The file to print parameters to. If None, the parameters are not printed,
            and this function simply acts as a test that all the parameters can be
            evaluated, and returns the dictionary of evaluated parameters.

        Returns
        -------
        evaluated_parameters : defaultdict
            The evaluated parameters, for further processing if needed

        Notes
        -----
        A C-rate of 1 C is the current required to fully discharge the battery in 1
        hour, 2 C is current to discharge the battery in 0.5 hours, etc
        """
        # Set list of attributes to ignore, for when we are evaluating parameters from
        # a class of parameters
        ignore = [
            "__name__",
            "__doc__",
            "__package__",
            "__loader__",
            "__spec__",
            "__file__",
            "__cached__",
            "__builtins__",
            "absolute_import",
            "division",
            "print_function",
            "unicode_literals",
            "pybamm",
            "_options",
            "constants",
            "np",
            "geo",
            "elec",
            "therm",
            "half_cell",
            "x",
            "r",
        ]

        # If 'parameters' is a class, extract the dict
        if not isinstance(parameters, dict):
            parameters_dict = {
                k: v for k, v in parameters.__dict__.items() if k not in ignore
            }
            for domain in ["n", "s", "p"]:
                domain_param = getattr(parameters, domain)
                parameters_dict.update(
                    {
                        f"{domain}.{k}": v
                        for k, v in domain_param.__dict__.items()
                        if k not in ignore
                    }
                )
            parameters = parameters_dict

        evaluated_parameters = defaultdict(list)
        # Calculate parameters for each C-rate
        for Crate in [1, 10]:
            # Update Crate
            capacity = self.get("Nominal cell capacity [A.h]")
            if capacity is not None:
                self.update(
                    {"Current function [A]": Crate * capacity},
                    check_already_exists=False,
                )

            # Turn to regular dictionary for faster KeyErrors
            self._dict_items = dict(self._dict_items)

            for name, symbol in parameters.items():
                if isinstance(symbol, pybamm.Symbol):
                    try:
                        proc_symbol = self.process_symbol(symbol)
                    except KeyError:
                        # skip parameters that don't have a value in that parameter set
                        proc_symbol = None
                    if not (
                        callable(proc_symbol)
                        or proc_symbol is None
                        or proc_symbol.has_symbol_of_classes(
                            (pybamm.Concatenation, pybamm.Broadcast)
                        )
                    ):
                        evaluated_parameters[name].append(proc_symbol.evaluate(t=0))

            # Turn back to FuzzyDict
            self._dict_items = pybamm.FuzzyDict(self._dict_items)

        # Calculate C-dependence of the parameters based on the difference between the
        # value at 1C and the value at C / 10
        for name, values in evaluated_parameters.items():
            if values[1] == 0 or abs(values[0] / values[1] - 1) < 1e-10:
                C_dependence = ""
            elif abs(values[0] / values[1] - 10) < 1e-10:
                C_dependence = " * Crate"
            elif abs(values[0] / values[1] - 0.1) < 1e-10:
                C_dependence = " / Crate"
            evaluated_parameters[name] = (values[0], C_dependence)
        # Print the evaluated_parameters dict to output_file
        if output_file:
            self.print_evaluated_parameters(evaluated_parameters, output_file)

        return evaluated_parameters

    def print_evaluated_parameters(self, evaluated_parameters, output_file):
        """
        Print a dictionary of evaluated parameters to an output file

        Parameters
        ----------
        evaluated_parameters : defaultdict
            The evaluated parameters, for further processing if needed
        output_file : string, optional
            The file to print parameters to. If None, the parameters are not printed,
            and this function simply acts as a test that all the parameters can be
            evaluated

        """
        # Get column width for pretty printing
        column_width = max(len(name) for name in evaluated_parameters.keys())
        s = "{{:>{}}}".format(column_width)
        with open(output_file, "w") as file:
            for name, (value, C_dependence) in sorted(evaluated_parameters.items()):
                if 0.001 < abs(value) < 1000:
                    file.write(
                        (s + " : {:10.4g}{!s}\n").format(name, value, C_dependence)
                    )
                else:
                    file.write(
                        (s + " : {:10.3E}{!s}\n").format(name, value, C_dependence)
                    )

    @staticmethod
    def find_parameter(path):
        """Look for parameter file in the different locations
        in PARAMETER_PATH
        """
        # Check for absolute path
        if os.path.isfile(path) and os.path.isabs(path):
            pybamm.logger.verbose(f"Using absolute path: '{path}'")
            return path
        for location in pybamm.PARAMETER_PATH:
            trial_path = os.path.join(location, path)
            if os.path.isfile(trial_path):
                pybamm.logger.verbose(f"Using path: '{location}' + '{path}'")
                return trial_path
        raise FileNotFoundError("Could not find parameter {}".format(path))
