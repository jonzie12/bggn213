Analysis of the 2019-nCov Outbreak
================

## Coronavirus

Here we analyze infection data for the 2019 novel Coronavirus COVID-19
(2019-nCoV) epidemic. The raw data is pulled from the Johns Hopkins
University Center for Systems Science and Engineering (JHU CCSE)
Coronavirus repository.

A CSV file is available here
<https://github.com/RamiKrispin/coronavirus-csv>

``` r
virus <- read.csv("coronavirus_dataset.csv")

tail(virus)
```

    ##      Province.State Country.Region     Lat     Long       date cases      type
    ## 2675         Shanxi Mainland China 37.5777 112.2922 2020-03-03     5 recovered
    ## 2676        Sichuan Mainland China 30.6171 102.7103 2020-03-03     8 recovered
    ## 2677        Tianjin Mainland China 39.3054 117.3230 2020-03-03    13 recovered
    ## 2678       Xinjiang Mainland China 41.1129  85.2401 2020-03-03     2 recovered
    ## 2679         Yunnan Mainland China 24.9740 101.4870 2020-03-03     1 recovered
    ## 2680       Zhejiang Mainland China 29.1832 120.0934 2020-03-03    24 recovered

Q1. How many total infected cases are there around the world?

``` r
A1<-sum(virus$cases)
print(A1)
```

    ## [1] 144233

Q2. How many deaths linked to infected cases have there been

``` r
A2<-sum(virus$cases[virus$type=="death"])
print(A2)
```

    ## [1] 3160

Q3. What is the overall death rate?

``` r
A3<-100*A2/A1
print(A3)
```

    ## [1] 2.190899

Q4. What is the death rate in China?

``` r
CasesinChina<- sum(virus$cases[virus$Country.Region=="Mainland China"])
DeathsinChina<- sum(virus$cases[virus$Country.Region=="Mainland China"&virus$type=="death"])
A4<-100*DeathsinChina/CasesinChina
A4
```

    ## [1] 2.256705

Q5. What is the death rate in Italy, Iran, and the US?

``` r
CasesinItaly<- sum(virus$cases[virus$Country.Region=="Italy"])
DeathsinItaly<- sum(virus$cases[virus$Country.Region=="Italy"&virus$type=="death"])
A5.Italy<-100*DeathsinItaly/CasesinItaly


CasesinIran<- sum(virus$cases[virus$Country.Region=="Iran"])
DeathsinIran<- sum(virus$cases[virus$Country.Region=="Iran"&virus$type=="death"])
A5.Iran<-100*DeathsinIran/CasesinIran


CasesinUS<- sum(virus$cases[virus$Country.Region=="US"])
DeathsinUS<- sum(virus$cases[virus$Country.Region=="US"&virus$type=="death"])
A5.US<-100*DeathsinUS/CasesinUS


print(A5.Iran)
```

    ## [1] 2.847633

``` r
print(A5.Italy)
```

    ## [1] 2.88216

``` r
print(A5.US)
```

    ## [1] 5.109489
