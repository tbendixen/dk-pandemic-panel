###################################################
### VELUX FOUNDATION COVID-19 RELIGIOSITY PANEL ###
### version 2.0, with four waves                ###
###################################################

################
### CODEBOOK ###
################

### FILE: "main_data_v2.csv"

## NB: YouGov turns "Don't know" into NA, but most items included a "Don't know" and/or
## "Prefer not to respond" option(s).

## id: unique id number

## wave: wave of data collection
# 1 14 - 21 May 2020
# 2 23 Sept. - 9. Oct. 2020
# 3 8 - 21 March 2021
# 4 3rd - 15th December 2021

## y (in wave 1)
# "Hvor enig eller uenig er du i udsagnet "Religion er meget vigtigt for mig"?"
# ("To what extent do you agree or disagree with the following: "Religion is very important to me"?")
# 1	Helt uenig                (strongly disagree)
# 2	Overvejende uenig         (somewhat disagree)
# 3	Hverken enig eller uenig  (neither agree nor disagree)
# 4	Overvejende enig          (somewhat agree)
# 5	Helt enig                 (strongly agree)

## y (in wave 2 and 3)
# "Hvor vigtigt er religion i dit liv?" ("How important is religion to you?")
# 1 Meget vigtigt         (very important)
# 2 Temmelig vigtigt      (somewhat important)
# 3 Ikke særlig vigtigt   (not really important)
# 4 Slet ikke vigtigt     (not at all important)

## health: 
# "Hvordan vil du beskrive dit nuværende helbred?" ("How would you describe your current health?")
# 1 Meget dårligt (very bad)
# 2 Dårligt       (bad)
# 3 Nogenlunde    (okay)
# 4 Godt          (good)
# 5 Meget godt    (very good)

## gender:
# 1 woman
# 2 man

## age: age in whole years

## edu: educational level
# 1 Grund-/folkeskole                                 (primary school)
# 2 Almen gymnasial uddannelse (studentereksamen/HF)  (high school/vocational training)
# 3 Erhvervsgymnasial uddannelse (HH/HTX/HHX)         (high school/vocational training)
# 4 Erhvervsfaglig uddannelse                         (high school/vocational training)
# 6 Kort videregående uddannelse under 3 år           (short cycle higher education (below three years))
# 7 Mellemlang videregående uddannelse 3-4 år         ((vocational) bachelors education (3-4 years))
# 8 Lang videregående uddannelse 5 år eller mere      (masters education (5 years or more))
# 6 Forskeruddannelse (f.eks. ph.d.)                  (doctoral degree) 

## household_income: annual household income (pre-tax)
# 1  Mindre end 100.000 kr. (less than 100.000 kr.)
# 2  100.000 - 199.999 kr.
# 3  200.000 - 299.999 kr.
# 4  300.000 - 399.999 kr.
# 5  400.000 - 499.999 kr.
# 6  500.000 - 599.999 kr.
# 7  600.000 - 699.999 kr.
# 8  700.000 - 799.999 kr.
# 9  800.000 - 899.999 kr.
# 10 900.000 - 999.999 kr.
# 11 1.000.000 kr. eller mere (1.000.000 kr. or more)

### FILE: "reliability_data.csv"

## id: unique id number

## ritual: 
# "Har du deltaget i en/flere religiøse handlinger (fx gudstjeneste, bøn, offer, prædikener, foredrag, undervisning osv.) siden jul?"
# ("Have you participated in one or more religious rituals or acts (e.g., church service, paryer, sacrifice, sermons, lectures, etc.) since Christmas ?")
# 1 Aldrig          (never)
# 2 Næsten aldrig   (almost never)
# 3 Ja, sjældent    (yes, rarely)
# 4 Ja, ofte        (yes, often)
# 5 Ved ikke        (don't know)

## church:
# "Hvor ofte plejer du at gå i kirke, moske, tempel, gurdwara, synagoge eller lignende?"
# ("How often do you usually attend church, mosque, temple, gurdwara, synagogue, or similar?")
# 1 Aldrig eller næsten aldrig  (never or almost never)
# 2 Flere gange om året         (several times a year)
# 3 Mindst én gang om måneden   (at least once a month)
# 4 Én gang om ugen             (once a week)
# 5 Mere end én gang om ugen    (more than once a week)
# 6 Hver dag                    (every day)

## belief:
# "Uanset om du går i kirke eller ej, vil du så mene du er..."
# ("Whether you attend church or not, would you say you are...")
# 1 Et troende menneske         ("a believer")
# 2 Et ikke troende menneske    ("a non-believer")
# 3 Overbevist ateist           ("confident atheist")
# 4 Ved ikke                    ("dont know")

## god:
# "Hvilket af disse udsagn kommer nærmest din tro?"
# ("Which of these statements are closest to your faith?")
# 1 Der er en personlig Gud                                               ("there's a personal God")
# 2 Der er en guddommelig kraft eller ånd, men ikke personificeret Gud    ("there's a divine force or spirit, but not a personified God")
# 3 Jeg ved ikke, hvad jeg skal tro                                       ("I don't know what to believe")
# 4 Jeg tror ikke der er nogen åndelig kraft eller personlig Gud          ("I don't believer there's a spiritual force or personal God")
# 5 Ved ikke                                                              ("don't know")

## y (in wave 1)
# "Hvor enig eller uenig er du i udsagnet "Religion er meget vigtigt for mig"?"
# ("To what extent do you agree or disagree with the following: "Religion is very important to me"?")
# 1	Helt uenig                (strongly disagree)
# 2	Overvejende uenig         (somewhat disagree)
# 3	Hverken enig eller uenig  (neither agree nor disagree)
# 4	Overvejende enig          (somewhat agree)
# 5	Helt enig                 (strongly agree)

### FILE: "exposure_data.csv"

## id: unique id number

## wave: wave of data collection
# 1 14 - 21 May 2020
# 2 23 Sept. - 9. Oct. 2020
# 3 8 - 21 March 2021

## y (in wave 1)
# "Hvor enig eller uenig er du i udsagnet "Religion er meget vigtigt for mig"?"
# ("To what extent do you agree or disagree with the following: "Religion is very important to me"?")
# 1	Helt uenig                (strongly disagree)
# 2	Overvejende uenig         (somewhat disagree)
# 3	Hverken enig eller uenig  (neither agree nor disagree)
# 4	Overvejende enig          (somewhat agree)
# 5	Helt enig                 (strongly agree)

## y (in wave 2 and 3)
# "Hvor vigtigt er religion i dit liv?" ("How important is religion to you?")
# 1 Meget vigtigt         (very important)
# 2 Temmelig vigtigt      (somewhat important)
# 3 Ikke særlig vigtigt   (not really important)
# 4 Slet ikke vigtigt     (not at all important)

### Pandemic exposure: 
# "Hvorledes har du haft Coronavirus inde på dit eget liv?"
# ("How has corona virus impacted your life/health?")

# Includes four sub items:

## pandemic_exo_1: pandemic exposure, self
# "Jeg har selv været syg" ("I've been ill [with corona virus] myself")
# 1                             Nej, men er ikke blevet testet (no, but have not been tested)
# 2                           Nej, og er blevet testet negativ (no, and have tested negative)
# 3 Ja, har testet positiv, men ikke særlig syg af Coronavirus (yes, have tested positive but was not very ill)
# 4                    Ja, har været meget syg med Coronavirus (yes, have tested positive and was very ill)
# 5                                                Ja, indlagt (yes, was hospitalized)

## pandemic_exo_2: pandemic exposure, household
# "Jeg har/har haft syge i min husstand" ("Someone in my household is/has been ill [with corona virus]")
# 1                             Nej, men er ikke blevet testet (no, but have not been tested)
# 2                           Nej, og er blevet testet negativ (no, and have tested negative)
# 3 Ja, har testet positiv, men ikke særlig syg af Coronavirus (yes, have tested positive but was not very ill)
# 4                    Ja, har været meget syg med Coronavirus (yes, have tested positive and was very ill)
# 5                                                Ja, indlagt (yes, was hospitalized)

## pandemic_exo_3: pandemic exposure, family
# "Jeg har/har haft syge i min family" ("Someone in my family is/has been ill [with corona virus]")
# 1                             Nej, men er ikke blevet testet (no, but have not been tested)
# 2                           Nej, og er blevet testet negativ (no, and have tested negative)
# 3 Ja, har testet positiv, men ikke særlig syg af Coronavirus (yes, have tested positive but was not very ill)
# 4                    Ja, har været meget syg med Coronavirus (yes, have tested positive and was very ill)
# 5                                                Ja, indlagt (yes, was hospitalized)

## pandemic_exo_4: pandemic exposure, relations
# "Jeg har/har haft syge i min bekendtskabskreds" ("Someone in my relations is/has been ill [with corona virus]")
# 1                             Nej, men er ikke blevet testet (no, but have not been tested)
# 2                           Nej, og er blevet testet negativ (no, and have tested negative)
# 3 Ja, har testet positiv, men ikke særlig syg af Coronavirus (yes, have tested positive but was not very ill)
# 4                    Ja, har været meget syg med Coronavirus (yes, have tested positive and was very ill)
# 5                                                Ja, indlagt (yes, was hospitalized)

