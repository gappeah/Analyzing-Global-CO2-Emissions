**Dataset**: CO2 Emissions by Country (Our World in Data)

- **Source**: [Our World in Data CO2 Dataset](https://ourworldindata.org/co2-and-greenhouse-gas-emissions)
- **Description**: Contains annual CO2 emissions (in million tonnes) by country, along with GDP, population, and per capita emissions from 1800 to recent years.
- **File**: Download the CSV file co2-data.csv.
https://github.com/owid/co2-data?tab=readme-ov-file

**Data Preparation**

- **Tasks**:
    - Load the dataset using Pandas.
    - Filter for years after 1990 and countries with non-null CO2 emissions.
    - Handle missing values (e.g., impute population using median or drop rows).
    - Create a new column for CO2 per capita using NumPy calculations.

**EDA**

- **Tasks**:
    - Use Pandas to compute summary statistics (mean, median CO2 emissions by continent).
    - Group by year to analyze global CO2 trends using NumPy for aggregations.
    - Identify top 5 emitting countries in 2020.

**Data Visualization**

- **Tasks**:
    - Plot a line chart of global CO2 emissions over time using Matplotlib.
    - Create a Seaborn boxplot of CO2 per capita by continent.
    - Use Seaborn to plot a scatterplot of CO2 vs. GDP with regression line.