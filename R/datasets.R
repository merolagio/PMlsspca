
#HS data=======================
#'  \strong{Holzinger-Swineford Student Ability data}
#'
#' These data are taken from the classic psychometrics dataset from Holzinger and Swineford (1939). The raw data are available from different sources, we used the R package \code{psychTools}, which also contain a detailed description of data and their historical usage. For comparability with previous analyses we selected 12 items and students from only from the Grant-White School (see also Ferrara, C., Martella, F., and Vichi, M. (2019). Probabilistic disjoint principal component analysis. Multivariate Behavioral Research, 54(1):47?61), 

#' @format 
#' A data frame with 145 rows and 12 items corresponding to four ability scales:
#' spatial (SPL)
#' \describe{
#' \item{visual}{Visual perception test, a nonlanguage multiple-choice test of spatial relations }
#' \item{cubes}{Cubes test, spatial relations}
#' \item{flags}{Lozenges test, a visual imagery test in two or three dimensions}
#' }
#' verbal (VBL)
#' \describe{
#' \item{paragraph}{Paragraph comprehension test, comprehension as measured by completion and multiple-choice questions}
#' \item{sentence}{Sentence completion test, a multiple-choice test in which “correct” answers reflect good judgment on the part of the subject}
#' \item{wordm}{Word meaning test, a multiple-choice vocabulary test}
#' }
#' speed (SPD)
#' \describe{
#' \item{addition}{Addition test, speed of adding pairs of one-digit numbers}
#' \item{counting}{Counting groups of dots test, 4–7 dots, arranged in random patterns to be counted by subject}
#' \item{straight}{Straight and curved capitals test, a series of capital letters to be distinguished between those composed of straight lines only and those containing curved lines}
#' }
#' mathematical (MTH)
#' \describe{
#' \item{deduct}{Deduction test, logical deduction test using the symbols and the letters}
#' \item{numeric}{Numerical puzzles test, a numerical deduction test, the object being to supply four numbers which will produce four given answers employing the operations of addition, multiplication, or division}
#' \item{series}{Series completion test, from a series of five numbers, the subject is supposed to deduce the rule of procedure from one number to the next and thus supply the sixth number in the series}
#' }
#' @details The data provided with this package are scaled to zero mean and unit variance.
#' @keywords datasets hs


#MSSCQ data=======================
#' Multidimensional Sexual Self-Concept Questionnaire (MSSCQ) data
#'
#' The data come from an interactive version of the MSSCQ hosted by OpenPsychometrics.
#' Sexual self-concept refers to a person's view of their own sexual behaviors and actions.
#' The MSSCQ was created by William E. Snell, Jr. (1998) for the general study of sexuality
#' and measures 20 scales.
#'
#' The dataset contains 17,685 responses to 100 statements, each rated on a 5-point Likert
#' scale: 1 = Not at all characteristic of me, 2 = Slightly characteristic of me,
#' 3 = Somewhat characteristic of me, 4 = Moderately characteristic of me,
#' 5 = Very characteristic of me. Items belonging to the same scale were placed together
#' in the questionnaire. The data in this package are ordered by scale membership. Scores
#' for the six reverse-scored items are reversed.
#'
#' @source OpenPsychometrics raw data: \url{https://openpsychometrics.org/_rawdata/};
#'   test page: \url{https://openpsychometrics.org/tests/MSSCQ.php}.
#'
#' @format A data frame with 17,685 rows (respondents) and 100 columns (items), coded on a
#'   1--5 Likert scale.
#'   
#' The table reports the item identifier, position in the questionnaire, the scale name,
#' a short scale code, and the full statement text. Items marked with `(R)` are reverse
#' items (as indicated in the statement text).
#'
#' @format A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{Item}{Item identifier (`I1`--`I100`).}
#'   \item{Position}{Question position in the questionnaire (`Q1`--`Q100`).}
#'   \item{Scale}{Scale name.}
#'   \item{Code}{Short scale code.}
#'   \item{Statement}{Item statement shown to respondents.}
#' }
#'
#' @section Items:
#' \tabular{lllll}{
#' \strong{Item} \tab \strong{Position} \tab \strong{Scale} \tab \strong{Code} \tab \strong{Statement} \cr
#' I1 \tab Q1 \tab sexual anxiety \tab SAN \tab I feel anxious when I think about the sexual aspects of my life. \cr
#' I2 \tab Q21 \tab sexual anxiety \tab SAN \tab I have the ability to take care of any sexual needs and desires that I may have. \cr
#' I3 \tab Q41 \tab sexual anxiety \tab SAN \tab I am very aware of my sexual feelings and needs. \cr
#' I4 \tab Q61 \tab sexual anxiety \tab SAN \tab I am motivated to avoid engaging in "risky" (i.e., unprotected) sexual behavior. \cr
#' I5 \tab Q81 \tab sexual anxiety \tab SAN \tab The sexual aspects of my life are determined mostly by chance happenings. \cr
#' I6 \tab Q2 \tab sexual self efficacy \tab SSE \tab I think about sex "all the time." \cr
#' I7 \tab Q22 \tab sexual self efficacy \tab SSE \tab I'm very assertive about the sexual aspects of my life. \cr
#' I8 \tab Q42 \tab sexual self efficacy \tab SSE \tab I expect that the sexual aspects of my life will be positive and rewarding in the future. \cr
#' I9 \tab Q62 \tab sexual self efficacy \tab SSE \tab I would be to blame, if the sexual aspects of my life were not going very well. \cr
#' I10 \tab Q82 \tab sexual self efficacy \tab SSE \tab I notice how others perceive and react to the sexual aspects of my life. \cr
#' I11 \tab Q3 \tab sexual consciousness \tab SC \tab I'm motivated to be sexually active. \cr
#' I12 \tab Q23 \tab sexual consciousness \tab SC \tab If I were to experience a sexual problem, I myself would in control of whether this improved. \cr
#' I13 \tab Q43 \tab sexual consciousness \tab SC \tab I derive a sense of self-pride from the way I handle my own sexual needs and desires. \cr
#' I14 \tab Q63 \tab sexual consciousness \tab SC \tab I am satisfied with the way my sexual needs are currently being met. \cr
#' I15 \tab Q83 \tab sexual consciousness \tab SC \tab My sexual behaviors are determined largely by other more powerful and influential people. \cr
#' I16 \tab Q4 \tab motivation to avoid risky sex \tab MTA \tab Not only would I be a good sexual partner, but it's quite important to me that I be a good sexual partner. \cr
#' I17 \tab Q24 \tab motivation to avoid risky sex \tab MTA \tab I am afraid of becoming sexual involved with another person. \cr
#' I18 \tab Q44 \tab motivation to avoid risky sex \tab MTA \tab If I am careful, then I will be able to prevent myself from having any sexual problems. \cr
#' I19 \tab Q64 \tab motivation to avoid risky sex \tab MTA \tab I am depressed about the sexual aspects of my life. \cr
#' I20 \tab Q84 \tab motivation to avoid risky sex \tab MTA \tab My sexuality is something that I am largely responsible for. \cr
#' I21 \tab Q5 \tab chance luck sexual control \tab CLS \tab I worry about the sexual aspects of my life. \cr
#' I22 \tab Q25 \tab chance luck sexual control \tab CLS \tab I am competent enough to make sure that my sexual needs are fulfilled. \cr
#' I23 \tab Q45 \tab chance luck sexual control \tab CLS \tab I am very aware of my sexual motivations and desires. \cr
#' I24 \tab Q65 \tab chance luck sexual control \tab CLS \tab I am motivated to keep myself from having any "risky" sexual behavior (e.g., exposure to sexual diseases). \cr
#' I25 \tab Q85 \tab chance luck sexual control \tab CLS \tab Most things that affect the sexual aspects of my life happen to me by accident. \cr
#' I26 \tab Q6 \tab sexual preoccupation \tab SP \tab I think about sex more than anything else. \cr
#' I27 \tab Q26 \tab sexual preoccupation \tab SP \tab I'm not very direct about voicing my sexual needs and preferences. (R) \cr
#' I28 \tab Q46 \tab sexual preoccupation \tab SP \tab I believe that in the future the sexual aspects of my life will be healthy and positive. \cr
#' I29 \tab Q66 \tab sexual preoccupation \tab SP \tab If the sexual aspects of my life were to go wrong, I would be the person to blame. \cr
#' I30 \tab Q86 \tab sexual preoccupation \tab SP \tab I'm concerned with how others evaluate my own sexual beliefs and behaviors. \cr
#' I31 \tab Q7 \tab sexual assertiveness \tab SAS \tab I'm motivated to devote time and effort to sex. \cr
#' I32 \tab Q27 \tab sexual assertiveness \tab SAS \tab If I were to experiences a sexual problem, my own behavior would determine whether I improved. \cr
#' I33 \tab Q47 \tab sexual assertiveness \tab SAS \tab I am proud of the way I deal with and handle my own sexual desires and needs. \cr
#' I34 \tab Q67 \tab sexual assertiveness \tab SAS \tab I am satisfied with the status of my own sexual fulfillment. \cr
#' I35 \tab Q87 \tab sexual assertiveness \tab SAS \tab My sexual behaviors are largely controlled by people other than myself (e.g., my partner, friends, family). \cr
#' I36 \tab Q8 \tab sexual optimism \tab SO \tab Not only would I be a skilled sexual partner, but it's very important to me that I be a skilled sexual partner. \cr
#' I37 \tab Q28 \tab sexual optimism \tab SO \tab I have a fear of sexual relationships. \cr
#' I38 \tab Q48 \tab sexual optimism \tab SO \tab I can pretty much prevent myself from developing sexual problems by taking good care of myself. \cr
#' I39 \tab Q68 \tab sexual optimism \tab SO \tab I am disappointed about the quality of my sex life. \cr
#' I40 \tab Q88 \tab sexual optimism \tab SO \tab The sexual aspects of my life are determined in large part by my own behavior. \cr
#' I41 \tab Q9 \tab sexual problem self blame \tab SPS \tab Thinking about the sexual aspects of my life often leaves me with an uneasy feeling. \cr
#' I42 \tab Q29 \tab sexual problem self blame \tab SPS \tab I have the skills and ability to ensure rewarding sexual behaviors for myself. \cr
#' I43 \tab Q49 \tab sexual problem self blame \tab SPS \tab I tend to think about my own sexual beliefs and attitudes. \cr
#' I44 \tab Q69 \tab sexual problem self blame \tab SPS \tab I want to avoid engaging in sex where I might be exposed to sexual diseases. \cr
#' I45 \tab Q89 \tab sexual problem self blame \tab SPS \tab Luck plays a big part in influencing the sexual aspects of my life. \cr
#' I46 \tab Q10 \tab sexual monitoring \tab SMN \tab I tend to be preoccupied with sex. \cr
#' I47 \tab Q30 \tab sexual monitoring \tab SMN \tab I am somewhat passive about expressing my own sexual desires. (R) \cr
#' I48 \tab Q50 \tab sexual monitoring \tab SMN \tab I do not expect to suffer any sexual problems or frustrations in the future. \cr
#' I49 \tab Q70 \tab sexual monitoring \tab SMN \tab If I were to develop a sexual disorder, then I would be to blame for not taking good care of myself. \cr
#' I50 \tab Q90 \tab sexual monitoring \tab SMN \tab I am quick to notice other people's reactions to the sexual aspects of my own life. \cr
#' I51 \tab Q11 \tab sexual motivation \tab SMT \tab I have a desire to be sexually active. \cr
#' I52 \tab Q31 \tab sexual motivation \tab SMT \tab If I were to become sexually maladjusted, I myself would be responsible for making myself better. \cr
#' I53 \tab Q51 \tab sexual motivation \tab SMT \tab I am pleased with how I handle my own sexual tendencies and behaviors. \cr
#' I54 \tab Q71 \tab sexual motivation \tab SMT \tab The sexual aspects of my life are personally gratifying to me. \cr
#' I55 \tab Q91 \tab sexual motivation \tab SMT \tab My sexual behavior is determined by the actions of powerful others (e.g., my partner, friends, family). \cr
#' I56 \tab Q12 \tab sexual problem management \tab SPM \tab Not only could I relate well to a sexual partner, but it's important to me that I be able to do so. \cr
#' I57 \tab Q32 \tab sexual problem management \tab SPM \tab I am fearful of engaging sexual activity. \cr
#' I58 \tab Q52 \tab sexual problem management \tab SPM \tab If just I look out for myself, then I will be able to avoid any sexual problems in the future. \cr
#' I59 \tab Q72 \tab sexual problem management \tab SPM \tab I feel discouraged about my sex life. \cr
#' I60 \tab Q92 \tab sexual problem management \tab SPM \tab I am in control of and am responsible for the sexual aspects of my life. \cr
#' I61 \tab Q13 \tab sexual esteem \tab SE \tab I worry about the sexual aspects of my life. \cr
#' I62 \tab Q33 \tab sexual esteem \tab SE \tab I am able to cope with and to handle my own sexual needs and wants. \cr
#' I63 \tab Q53 \tab sexual esteem \tab SE \tab I'm very alert to changes in my sexual thoughts, feelings, and desires. \cr
#' I64 \tab Q73 \tab sexual esteem \tab SE \tab I really want to prevent myself from being exposed to sexual diseases. \cr
#' I65 \tab Q93 \tab sexual esteem \tab SE \tab The sexual aspects of my life are largely a matter of (good or bad) fortune. \cr
#' I66 \tab Q14 \tab sexual satisfaction \tab SS \tab I'm constantly thinking about having sex. \cr
#' I67 \tab Q34 \tab sexual satisfaction \tab SS \tab I do not hesitate to ask for what I want in a sexual relationship. \cr
#' I68 \tab Q54 \tab sexual satisfaction \tab SS \tab I will probably experience some sexual problems in the future. (R) \cr
#' I69 \tab Q74 \tab sexual satisfaction \tab SS \tab If I were to develop a sexual problem, then it would be my own fault for letting it happen. \cr
#' I70 \tab Q94 \tab sexual satisfaction \tab SS \tab I'm concerned about how the sexual aspects of my life appear to others. \cr
#' I71 \tab Q15 \tab power other sexual control \tab POS \tab It's important to me that I involve myself in sexual activity. \cr
#' I72 \tab Q35 \tab power other sexual control \tab POS \tab If I developed any sexual problems, my recovery would depend in large part on what I myself would do. \cr
#' I73 \tab Q55 \tab power other sexual control \tab POS \tab I have positive feelings about the way I approach my own sexual needs and desires. \cr
#' I74 \tab Q75 \tab power other sexual control \tab POS \tab The sexual aspects of my life are satisfactory, compared to most people's. \cr
#' I75 \tab Q95 \tab power other sexual control \tab POS \tab In order to be sexually active, I have to conform to other more powerful individuals. \cr
#' I76 \tab Q16 \tab sexual self schemata \tab SSS \tab I am able to "connect" well with a sexual partner, and it's important to me that I am able to do so. \cr
#' I77 \tab Q36 \tab sexual self schemata \tab SSS \tab I don't have much fear about engaging in sex. (R) \cr
#' I78 \tab Q56 \tab sexual self schemata \tab SSS \tab I will be able to avoid any sexual problems, if I just take good care of myself. \cr
#' I79 \tab Q76 \tab sexual self schemata \tab SSS \tab I feel unhappy about my sexual experiences. \cr
#' I80 \tab Q96 \tab sexual self schemata \tab SSS \tab The main thing which affects the sexual aspects of my life is what I myself do. \cr
#' I81 \tab Q17 \tab fear of sex \tab FOS \tab I feel nervous when I think abut the sexual aspects of my life. \cr
#' I82 \tab Q37 \tab fear of sex \tab FOS \tab I have the capability to take care of my own sexual needs and desires. \cr
#' I83 \tab Q57 \tab fear of sex \tab FOS \tab I am very aware of the sexual aspects of myself (e.g. habits, thoughts, beliefs). \cr
#' I84 \tab Q77 \tab fear of sex \tab FOS \tab I am really motivated to avoid any sexual activity that might expose me to sexual diseases. \cr
#' I85 \tab Q97 \tab fear of sex \tab FOS \tab The sexual aspects of my life are a matter of fate (destiny). \cr
#' I86 \tab Q18 \tab sexual problem prevention \tab SPP \tab I think about sex the majority of the time. \cr
#' I87 \tab Q38 \tab sexual problem prevention \tab SPP \tab When it comes to sex, I usually ask for what I want. \cr
#' I88 \tab Q58 \tab sexual problem prevention \tab SPP \tab I anticipate that in the future the sexual aspects of my life will be frustrating. (R) \cr
#' I89 \tab Q78 \tab sexual problem prevention \tab SPP \tab If something went wrong with my own sexuality, then it would be my own fault. \cr
#' I90 \tab Q98 \tab sexual problem prevention \tab SPP \tab I'm aware of the public impression created by my own sexual behaviors and attitudes. \cr
#' I91 \tab Q19 \tab sexual depression \tab SD \tab I strive to keep myself sexually active. \cr
#' I92 \tab Q39 \tab sexual depression \tab SD \tab If I developed a sexual disorder, my recovery would depend on how I myself dealt with the problem. \cr
#' I93 \tab Q59 \tab sexual depression \tab SD \tab I feel good about the way I express my own sexual needs and desires. \cr
#' I94 \tab Q79 \tab sexual depression \tab SD \tab I am satisfied with the sexual aspects of my life. \cr
#' I95 \tab Q99 \tab sexual depression \tab SD \tab My sexual behavior is mostly determined by people who have influence and control over me. \cr
#' I96 \tab Q20 \tab internal sexual control \tab ISC \tab Not only am I be capable of relating to a sexual partner, but it's important to me that I relate very well. \cr
#' I97 \tab Q40 \tab internal sexual control \tab ISC \tab I'm not afraid of becoming sexually active. (R) \cr
#' I98 \tab Q60 \tab internal sexual control \tab ISC \tab If I just pay careful attention, I'll be able to prevent myself from having any sexual problems. \cr
#' I99 \tab Q80 \tab internal sexual control \tab ISC \tab I feel sad when I think about my sexual experiences. \cr
#' I100 \tab Q100 \tab internal sexual control \tab ISC \tab My sexuality is something that I myself am in charge of.
#' }
#'
#' @keywords datasets msscq

# Crime data===============
#' The \code{cr} dataset combines socio-economic data from the 1990 US Census, law enforcement data from the 1990 US LEMAS survey, and crime data from the 1995 FBI UCR.
#'
#' In this package, the variable order was rearranged and 25 variables with missing values were removed. The variables in \code{cr} are scaled to mean zero and unit variance. The lookup data frame \code{crime_var_lookup_df} is also included.
#'
#' @source \url{http://archive.ics.uci.edu/ml/machine-learning-databases/communities/}
#'
#' @references
#' Redmond, M. A. and Baveja, A. (2002). A Data-Driven Software Tool for Enabling Cooperative
#' Information Sharing Among Police Departments. \emph{European Journal of Operational Research},
#' 141, 660--678.
#'
#' @format A numeric data frame with 1,994 observations and 99 variables scaled to mean zero and unit variance). See
#'   \code{crime_var_lookup_df} for variable definitions.
#'
#' @section Variables:
#' \tabular{ll}{
#' \strong{Code} \tab \strong{description} \cr
#' population \tab Population of the community. \cr
#' numbUrban \tab Number of residents living in areas classified as urban. \cr
#' pctUrban \tab Percent of residents living in areas classified as urban. \cr
#' PctBornSameState \tab Percent of residents born in the same state where they currently live. \cr
#' PctSameHouse85 \tab Percent of residents living in the same house as in 1985. \cr
#' PctSameCity85 \tab Percent of residents living in the same city as in 1985. \cr
#' PctSameState85 \tab Percent of residents living in the same state as in 1985. \cr
#' LandArea \tab Land area (square miles). \cr
#' PopDens \tab Population density (persons per square mile). \cr
#' PctUsePubTrans \tab Percent of residents using public transit to commute. \cr
#' LemasPctOfficDrugUn \tab Percent of officers assigned to drug units. \cr
#' racepctblack \tab Percent of residents who are African American. \cr
#' racePctWhite \tab Percent of residents who are White. \cr
#' racePctAsian \tab Percent of residents who are of Asian heritage. \cr
#' racePctHisp \tab Percent of residents who are of Hispanic heritage. \cr
#' agePct12t21 \tab Percent of residents aged 12--21. \cr
#' agePct12t29 \tab Percent of residents aged 12--29. \cr
#' agePct16t24 \tab Percent of residents aged 16--24. \cr
#' agePct65up \tab Percent of residents aged 65 and older. \cr
#' NumImmig \tab Total number of residents known to be foreign born. \cr
#' PctImmigRecent \tab Percent of immigrants who immigrated within the last 3 years. \cr
#' PctImmigRec5 \tab Percent of immigrants who immigrated within the last 5 years. \cr
#' PctImmigRec8 \tab Percent of immigrants who immigrated within the last 8 years. \cr
#' PctImmigRec10 \tab Percent of immigrants who immigrated within the last 10 years. \cr
#' PctRecentImmig \tab Percent of residents who immigrated within the last 3 years. \cr
#' PctRecImmig5 \tab Percent of residents who immigrated within the last 5 years. \cr
#' PctRecImmig8 \tab Percent of residents who immigrated within the last 8 years. \cr
#' PctRecImmig10 \tab Percent of residents who immigrated within the last 10 years. \cr
#' PctSpeakEnglOnly \tab Percent of residents who speak only English. \cr
#' PctNotSpeakEnglWell \tab Percent of residents who do not speak English well. \cr
#' PctForeignBorn \tab Percent of residents who are foreign born. \cr
#' MalePctDivorce \tab Percent of males who are divorced. \cr
#' MalePctNevMarr \tab Percent of males who have never married. \cr
#' FemalePctDiv \tab Percent of females who are divorced. \cr
#' TotalPctDiv \tab Percent of the population who are divorced. \cr
#' PersPerFam \tab Mean number of persons per family. \cr
#' PctFam2Par \tab Percent of families with children headed by two parents. \cr
#' PctKids2Par \tab Percent of children in family housing living with two parents. \cr
#' PctYoungKids2Par \tab Percent of children aged 4 and under living with two parents. \cr
#' PctTeen2Par \tab Percent of children aged 12--17 living with two parents. \cr
#' PctWorkMomYoungKids \tab Percent of mothers of children aged 6 and under in the labor force. \cr
#' PctWorkMom \tab Percent of mothers of children under 18 in the labor force. \cr
#' NumIlleg \tab Number of children born to never-married mothers. \cr
#' PctIlleg \tab Percent of children born to never-married mothers. \cr
#' householdsize \tab Mean number of persons per household. \cr
#' PctLargHouseFam \tab Percent of family households that are large (6 or more members). \cr
#' PctLargHouseOccup \tab Percent of occupied households that are large (6 or more persons). \cr
#' PersPerOccupHous \tab Mean persons per occupied household. \cr
#' PersPerOwnOccHous \tab Mean persons per owner-occupied household. \cr
#' PersPerRentOccHous \tab Mean persons per renter-occupied household. \cr
#' PctPersOwnOccup \tab Percent of persons living in owner-occupied households. \cr
#' PctPersDenseHous \tab Percent of persons in dense housing (more than 1 person per room). \cr
#' PctHousLess3BR \tab Percent of housing units with fewer than 3 bedrooms. \cr
#' MedNumBR \tab Median number of bedrooms. \cr
#' HousVacant \tab Number of vacant housing units. \cr
#' PctHousOccup \tab Percent of housing units that are occupied. \cr
#' PctHousOwnOcc \tab Percent of households that are owner occupied. \cr
#' PctVacantBoarded \tab Percent of vacant housing units that are boarded up. \cr
#' PctVacMore6Mos \tab Percent of vacant housing units vacant for more than 6 months. \cr
#' MedYrHousBuilt \tab Median year housing units were built. \cr
#' PctHousNoPhone \tab Percent of occupied housing units without a telephone. \cr
#' PctWOFullPlumb \tab Percent of housing units without complete plumbing facilities. \cr
#' OwnOccLowQuart \tab Owner-occupied housing value (lower quartile). \cr
#' OwnOccMedVal \tab Owner-occupied housing value (median). \cr
#' OwnOccHiQuart \tab Owner-occupied housing value (upper quartile). \cr
#' RentLowQ \tab Rent (lower quartile). \cr
#' RentMedian \tab Rent (median). \cr
#' RentHighQ \tab Rent (upper quartile). \cr
#' MedRent \tab Median gross rent (includes utilities). \cr
#' MedRentPctHousInc \tab Median gross rent as a percent of household income. \cr
#' MedOwnCostPctInc \tab Median owner cost as a percent of household income (with mortgage). \cr
#' MedOwnCostPctIncNoMtg \tab Median owner cost as a percent of household income (without mortgage). \cr
#' NumInShelters \tab Number of people in homeless shelters. \cr
#' NumStreet \tab Number of homeless people counted on the street. \cr
#' medIncome \tab Median household income. \cr
#' medFamInc \tab Median family income. \cr
#' perCapInc \tab Per capita income. \cr
#' whitePerCap \tab Per capita income for White residents. \cr
#' blackPerCap \tab Per capita income for African American residents. \cr
#' indianPerCap \tab Per capita income for Native American residents. \cr
#' AsianPerCap \tab Per capita income for residents of Asian heritage. \cr
#' HispPerCap \tab Per capita income for residents of Hispanic heritage. \cr
#' NumUnderPov \tab Number of residents below the poverty level. \cr
#' PctPopUnderPov \tab Percent of residents below the poverty level. \cr
#' PctLess9thGrade \tab Percent of residents aged 25+ with less than a 9th-grade education. \cr
#' PctNotHSGrad \tab Percent of residents aged 25+ who are not high school graduates. \cr
#' PctBSorMore \tab Percent of residents aged 25+ with a bachelor's degree or higher. \cr
#' PctUnemployed \tab Percent of residents aged 16+ in the labor force who are unemployed. \cr
#' PctEmploy \tab Percent of residents aged 16+ who are employed. \cr
#' PctEmplManu \tab Percent of residents aged 16+ employed in manufacturing. \cr
#' PctEmplProfServ \tab Percent of residents aged 16+ employed in professional services. \cr
#' PctOccupManu \tab Percent of residents aged 16+ in manufacturing occupations. \cr
#' PctOccupMgmtProf \tab Percent of residents aged 16+ in management or professional occupations. \cr
#' pctWWage \tab Percent of households with wage or salary income (1989). \cr
#' pctWFarmSelf \tab Percent of households with farm or self-employment income (1989). \cr
#' pctWInvInc \tab Percent of households with investment or rental income (1989). \cr
#' pctWSocSec \tab Percent of households with Social Security income (1989). \cr
#' pctWPubAsst \tab Percent of households with public assistance income (1989). \cr
#' pctWRetire \tab Percent of households with retirement income (1989). \cr
#' }
#' @keywords dataset cr




