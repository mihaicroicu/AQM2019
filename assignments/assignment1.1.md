## ASSIGNMENT 1.1 ## 

**Please do all of the EASY AND MEDIUM tasks**

**Take your pick at at least 1 (ONE) of the HARD tasks**

In order to pass, you do not need to successfully complete the HARD tasks, but at least have tried one of them.

Two commonly used datasets in conflict research are the UCDP Non-State Conflict dataset, which provides a list of armed conflicts where *states were not involved* on a country-year level and the UCDP Geo-Referenced Event Dataset (GED) which provides details on individual battles (down to day and village level). GED provides battle events for all three conflict types, *state-based*, *non-state* and *one-sided*. 

I've prepared two simplified and cleaned-up versions of the two datasets here:

GED: https://raw.githubusercontent.com/mihaicroicu/AQM2019/master/assignments/data/ged_assign1.csv
Non-State: https://raw.githubusercontent.com/mihaicroicu/AQM2019/master/assignments/data/nonstate_assign1.csv

Please use these files *and only these files* for the assignment.

A codebook describing what is in these files is available here: https://github.com/mihaicroicu/AQM2019/blob/master/assignments/codebook.md


Tasks:

Explore the two datasets and write scripts doing the following:

1. `EASY` `code` Calculate a yearly count of battle events and a yearly sum of total fatalites for the *state based* category. State-based conflict is conflict having a `type_of_violence == 1`. Provide the code doing this.

2. `EASY` `code` How many distinct state-based conflicts (not events, but conflicts) are **NOT** active in each year in the GED dataset? Provide the code doing this.

A common problem researchers and anlysts have is that non-state conflict is sub-divided into three sub-types called organization types (communal, organized and supportive). This data is not available in `GED` but is available in `Non-State`. You want to know the behavior of those conflicts based around a *communal* organization type (this is called "communal conflict" or "tribal conflict in literature).

Thus:

3. `MEDIUM` `code` Write the code that produces a dataset containing **sums of fatalities** and **counts of battle events** by **conflict-year** and **only covers communally organized non-state conflict**. For your research you need only **active** conflict-years included in your dataset. 

4. `MEDIUM` `write` Write a short (1-3 paragraphs) summary describing *the choices you made* and the *plan* you did to solve problem **3.**

Another common problem researchers have is needing monthly data. You have this problem for your research too!

Thus:

6. `MEDIUM` `write` Make a dataset where you summarize `GED` computing **sums of fatalities** and **counts of battle events** by country-month only for **Africa**, only for **2016** and only for **state-based conflict**.

**CHOOSE ONE BELOW**

7a. `HARD` `write+code` Your reviewer asks you to compute the **1-month lag** of **counts of battle events** and **sums of fatalities** for the dataset you computed at task **6.**. Can you do it using the dataset you computed at task **6.**? Why? Why not? Think about what happens to the lags when a month has no events and no fatalities in that month, but the month before and after does have events. Write 1-2 paragraphs. Does this dataset help you

7b. `HARD` `write+code` A conflict consists of two actors fighting each-other. 
