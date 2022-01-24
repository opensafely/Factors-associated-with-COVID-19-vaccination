# Factors-associated-with-COVID-19-vaccination

This is the code and configuration for the "Factors associated with" analysis wihtin the initial pre-print paper available on the [OpenSAFELY website here](https://opensafely.org/research/2021/covid-vaccine-coverage/) and on [MedRxiv here](https://www.medrxiv.org/content/10.1101/2021.01.25.21250356v2). 

## How to use
- If you are interested in how we defined our variables, take a look at the [study definition](analysis/study_definition_delivery.py); this is written in `python`, but non-programmers should be able to understand what is going on there.
- If you are interested in how we defined our code lists, look in the [codelists folder](./codelists/). All codelists are available online at [OpenCodelists](https://codelists.opensafely.org/) for inspection and re-use by anyone 
- Developers and epidemiologists interested in the framework should review [the OpenSAFELY documentation](https://docs.opensafely.org)


# About the OpenSAFELY framework

The OpenSAFELY framework is a secure analytics platform for
electronic health records research in the NHS.

Instead of requesting access for slices of patient data and
transporting them elsewhere for analysis, the framework supports
developing analytics against dummy data, and then running against the
real data *within the same infrastructure that the data is stored*.
Read more at [OpenSAFELY.org](https://opensafely.org).
