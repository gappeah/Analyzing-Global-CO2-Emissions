# YouTube Market Data Analysis

## Project Overview
This project is a comprehensive data analysis workflow that demonstrates how to transform raw data from Excel into actionable insights using SQL and Power BI. The dataset used contains information about top-performing YouTube channels in the UK, with data sourced from a publicly available platform.

The workflow consists of data extraction, cleaning, transformation, and visualization, designed to simulate a real-world marketing analysis scenario for a client aiming to identify YouTube influencers for marketing campaigns.

## Table of Contents
1. [Project Architecture](#project-architecture)
2. [Dataset Description](#dataset-description)
3. [Steps Involved](#steps-involved)
4. [Data Cleaning](#data-cleaning)
5. [Data Visualization](#data-visualization)
6. [Recommendations](#recommendations)
7. [How to Run](#how-to-run)
8. [Contributing](#contributing)
9. [License](#license)

## Project Architecture
The project involves three key technologies:
- **Excel**: The raw data is sourced from an Excel file containing YouTube channel information.
- **SQL**: SQL Server is used to clean, transform, and perform quality checks on the data.
- **Power BI**: A dashboard is created to visualize insights, including top YouTubers by subscriber count, views, and video uploads.

## Dataset Description
The dataset contains the following fields:
- Channel name
- Subscriber count
- Views
- Videos uploaded

The dataset was sourced from a data platform and extracted via a Python script using the YouTube API.

## Steps Involved

1. **Data Extraction**: 
   - The initial data is downloaded from an external source in Excel format.
   - A Python script is used to pull supplementary data from the YouTube API.

2. **Data Cleaning**:
   - Unnecessary columns are removed.
   - Data is cleaned to ensure consistency and accuracy.

3. **Data Transformation**:
   - Data types for each column are validated.
   - SQL queries are used to transform the dataset into the required format for analysis.

4. **Data Quality Checks**:
   - Row and column counts are verified.
   - Data types are validated.
   - Duplicate records are checked.

5. **Data Visualization**:
   - The cleaned data is imported into Power BI.
   - A dashboard is created with visualizations such as tables, tree maps, and bar charts to show top YouTubers by subscriber count, video uploads, and views.
