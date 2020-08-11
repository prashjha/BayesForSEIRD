# Steps to get csv data from excel file

- Remove first two colums

- Remove column at the very bottom

- Remove extra spacing between last and second last row

- Inster 0, 1, 2, ...,N row and column at the beginning. This will act as key for pandas dataframe.

- Change DeWitt to De Witt

- remove total population column for all three datas; total infected, active infected and deceased


## For 11 August 2020 data

- deceased data starts from 7 March 2020  
- infected data starts from 4 March 2020 and does not have entries for 7,8,14 March
- active data starts from 7 April
