# File information

- data files are in `data/` directory

- mesh files are in `data/mesh/` directory

Information about various files is provided below

## directory: mesh

- for coarse mesh use `mesh/mesh_10h` (924 vertices)

- for fine mesh use `mesh/mesh_4h` (4460 vertices) or `mesh_5h` (2969 vertices)

- mesh_2h has 15695 vertices

- mesh_h has 59396 vertices

- contains list of elements indices which belong to different districts in file `elements_in_district_mesh_10h.txt`

## directory: geography

- contains geography data which is fixed

### county_population.txt

- Total population in various counties

- Line 1 corresponds to population in county 1, 2 corresponds to population in county 2, and so on

- Name of counties for given county id can be found in file `county_names.txt`

### county_names.txt

- Names of counties

### county_name_and_population.txt

- Name and population of county side by side

### county_geom_details.txt

- county's centroid, area, and length

- for county i, the vector of data in line i+1 correspond to centroid's x coordinate, y coordinate, county's area and length

- load this file using: `np.loadtxt('county_geom_details.txt')`


### district_names.txt


- Name of districts

- For given district number i, the name is in the line i+1 of file

- We follow consistent naming and index convention in all files


### district_counties.txt

- Ids of counties belonging to district

- For district i, the ids of counties are in line i+1 of file

- Use following function to read the file

```py
N_districts = 25
dist_county = [0 for i in range(N_districts)]

fo = open(fdir + 'district_counties.txt', 'r')    
lcounty = fo.readline().split()
cnt = 1
while cnt <= N_districts:

    counties = []
    for s in lcounty:
        counties.append(int(s))
        
    dist_county[cnt-1] = counties
    
    lcounty = fo.readline().split()
    cnt += 1
```

## directory: data/covid_25June2020

- covid data upto 25 June 2020

### cases_data.txt

- day 0 is `6 March 2020`

- this file gives date for given day starting from day 0 upto last day of data

### deceased_county.txt

- deceased cases countywise

- for day i, the data is in line i+1 of file

- for day i, the data is the vector of cases where jth element in the vector correspond to case in county j

- load this file simply by `np.loadtxt('deceased_county.txt')`

### deceased_district.txt

- deceased cases districtwise

- for day i, the data is in line i+1 of file

- for day i, the data is the vector of cases where jth element in the vector correspond to case in district j

- load this file using `np.loadtxt('deceased_district.txt')`

### deceased_state.txt

- total deceased cases in state of Texas

- each row correspond to total cases at corresponding day

- load using: `np.loadtxt()`

### deceased_transition_day_county.txt

- day when county deceased cases transition from 0 to nonzero number

- those with -1 value mean that county cases have never exceeded 0

- load using `np.loadtxt()`

### infected_total_county.txt, infected_total_district.txt, infected_total_state.txt, infected_total_transition_day_county.txt

- same as files `deceased_district.txt` etc but now for total infected cases

- this is total infected cases which consists of active infected cases + recovered cases + deceased cases

- total infected cases and deceased cases are more accurate whereas recovered and active infected cases are estimates

### infected_active_county.txt, infected_active_district.txt, infected_active_state.txt, infected_active_transition_day_county.txt

- same as files `deceased_district.txt` etc but now for active infected cases

- this is active infected cases

### recovered_county.txt, recovered_district.txt, recovered_state.txt, recovered_transition_day_county.txt

- same as files `deceased_district.txt` etc but now for recovered cases

