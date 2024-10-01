**********************************************
**********************************************

** EXPLORATION AND STATISTICAL ANALYSIS ** 
** Paper: Booster association with neutralization and anti-S antibody titers in SARS-CoV-2 vaccinated adults in Mexico: a comparison across five types of vaccines.


** Mtro. Nabetse Baruc Blas Miranda
** Dra. Martha Carnalla
** Dr. Tonatiuh Barrientos
** National Institute of Public Health (INSP) of Mexico


** January, 2024

**********************************************
** Prepare database
**********************************************

clear
use "C:\Users\Nabetse Blas\OneDrive - Instituto Nacional de Salud Pública\Manuscripts\Paper GMTs\Paper GMTs FINAL\Comments Journal Med Virology\Base de datos paper GMTs\integrantes_ensanut2022_w copy.dta"
 

merge 1:1 upm est_sel FOLIO_I FOLIO_INT using "C:\Users\Nabetse Blas\OneDrive - Instituto Nacional de Salud Pública\Manuscripts\Paper GMTs\Paper GMTs FINAL\Comments Journal Med Virology\Base de datos paper GMTs\ensasangre22.sangre_rec_protn_30032023_5w copy.dta"

drop if _merge ==1
drop _merge

merge m:1 FOLIO_I using "C:\Users\Nabetse Blas\OneDrive - Instituto Nacional de Salud Pública\Manuscripts\Paper GMTs\Paper GMTs FINAL\Comments Journal Med Virology\Base de datos paper GMTs\ENSANUT2022_ICB copy.dta"

drop if _merge ==2
drop _merge

** subsample
drop if Interpretacion_S_PF==" "
** ((6,161  observations deleted)

**********************************************
** Cleaning data base
**********************************************

** <18 years old not elegible for this study
drop if h0303 <18
** (2,594 observations deleted)

** Participants unvaccinated 
drop if h1602==2
** (585 observations deleted)

** Participants who received J&J, Moderna or not remenber in 1st dose (Primary Scheme)
drop if h16051==6 | h16051==7 | h16051==9
** (313 observations deleted)

**********************************************

** Complete scheme // Complete scheme was defined as having one dose in CanSino and two doses were recommended in Pfizer, Sputnik, Sinovac and AztraZeneca, incomplete scheme otherwise.
keep if h16051==h16052 | h16051==5
** (594  observations deleted)

** Exclude fourth dose for all vaccines // Our objective includes only one booster //
drop if inrange(h16054, 1,9) 
* (282 observations deleted)

** Exclude third dose for CanSino 
drop if h16051==5 & inrange(h16053, 1,9)
* (73 observations deleted)

** Exclude not knowing or not remember in 2nd dose
drop if h16052==9
** (4 observations deleted)

** Exclude not knowing or not remember in 3rd dose
drop if h16053==9
*(16 observations deleted)

** Positivity anti-S
drop if Interpretacion_S_PF=="NEGATIVO"
** (3 observations deleted)


**********************************************
** Prepare anti-S variable
**********************************************

** Change relative units to international units of prot S
gen Resultado_S_PF_nueva= Resultado_S_PF*1.0149589261779

** Variable logarithm anti-S
gen var_log= ln(Resultado_S_PF_nueva)



**********************************************
** Create variables
**********************************************

** protN 
gen protN=.
replace protN=1 if Interpretacion_N_PF=="POSITIVO"
replace protN=0 if Interpretacion_N_PF=="NEGATIVO"

** sex
tab h0302
rename h0302 sexo

** age category in decades
gen edadcat=h0303
recode edadcat 18/29=1 30/39=2 40/49=3 50/59=4 60/100=5
label define edadcat 1"18-29" 2"30-39" 3"40-49" 4"50-59" 5"60+"
label values edadcat edadcat
tab edadcat

** education
gen escolaridad=h0317a
recode escolaridad 0/1=2 6=2 7=3 5=4 8=4 9/12=5
replace escolaridad=2 if escolaridad==. & inrange(h0303,1,2)
label define escolaridad 2 "Primaria o menos" 3"Secundaria" 4"Preparatoria" 5"Licenciatura o más" 
label values escolaridad escolaridad

** region // previously defined
tab region

** urbanization // previously defined
tab estrato

** socioeconomic status // previously defined
tab nseF

** security healthcare
global seguro   = "H0310A H0310B H0310C"
destring $seguro, replace
gen seguridad=0
foreach var in $seguro {
replace seguridad=1 if (inrange(`var',1,6) |`var'==9)
}

** employment status // only in participants >15 years 
gen ocupacion=h0324
recode ocupacion 1/7=3
replace ocupacion=0 if h0323==1
replace ocupacion=0 if h0323==2
replace ocupacion=1 if h0323==3
replace ocupacion=2 if h0323==4
replace ocupacion=0 if h0323==5
replace ocupacion=0 if h0323==6

sort h0324e
br h0324e if h0324==7

replace ocupacion=0 if h0324e=="AMA DE CASA"
replace ocupacion=0 if h0324e=="AYUDA A CUIDAR A SUS HIJAS"
replace ocupacion=0 if h0324e=="CUIDA A SUS POLLOS"
replace ocupacion=0 if h0324e=="DESCANSA"
replace ocupacion=0 if h0324e=="DESEMPLEADO"
replace ocupacion=0 if h0324e=="TRABAJO COSECHANDO MAIZ PARA SU PROPIO CONSUMO"

replace ocupacion=. if h0303<15
replace ocupacion=4 if ocupacion==3 & seguridad==1
label define ocupacion 0 "Desempleado" 1 "Estudiante" 2 "Jubilado/pensionado" 3 "Empleado_informal" 4"Empleado_formal"
label values ocupacion ocupacion

**********************************************
** indicator variable
**********************************************

**Type plataform
*Type of vaccines // code

** Messenger RNA
* Pfizer // 1
* Moderna // 6

** Viral vector
* Sputnik // 2
* Aztra Zeneca // 4
* CanSino // 5
* Johnson & Johnson // 7

** Inactivated virus
* Sinovac // 3

**********************************************
**********************************************

** Pfizer
gen pfizer=0 if h16051==1 & h16053==.
replace pfizer=1 if h16051==1 & h16053==1
replace pfizer=2 if h16051==1 & h16053==4
replace pfizer=3 if h16051==1 & h16053==2
replace pfizer=4 if h16051==1 & h16053==6
replace pfizer=5 if h16051==1 & h16053==5
replace pfizer=6 if h16051==1 & h16053==3
replace pfizer=7 if h16051==1 & h16053==7

label define pfizer 0"not booster" 1"pfizer" 2"aztra" 3"sputnik" 4"moderna" 5"Cansino" 6"sinovac" 7"JJ" , modify
label values pfizer pfizer

gen pfizermodelo=pfizer 
recode pfizermodelo (0=4) (4=2) (2=3) (5=3) (7=3) (6=.) 
***** Note: excluded Sinovac vaccine (not estimable)
label define booster  1"homologocis" 2"homologotrans" 3"heterologo" 4"sin booster" 5"herologoRNA", modify
***** Recod for colapse homologous
recode pfizermodelo (2=1)
label values pfizermodelo booster 
tab pfizermodelo pfizer


** Sputnik
gen sputnik=0 if h16051==2 & h16053==.
replace sputnik=1 if h16051==2 & h16053==2
replace sputnik=2 if h16051==2 & h16053==4
replace sputnik=3 if h16051==2 & h16053==1
replace sputnik=4 if h16051==2 & h16053==5
replace sputnik=5 if h16051==2 & h16053==6
replace sputnik=6 if h16051==2 & h16053==3

label define sputnik 0"not booster" 1"sputnik" 2"aztra" 3"pfizer" 4"Cansino" 5"moderna" 6"sinovac", modify
label values sputnik sputnik

** prepare varibale for adjusted model 
gen sputnikmodelo=sputnik 
recode sputnikmodelo (0=4) (4=2) (3=.) (5=.) (6=.)
***** Recod for colapse homologous
recode sputnikmodelo (2=1)
label values sputnikmodelo booster 
tab sputnikmodelo sputnik
***** Note. pfizer, cansino and sinovac (not estimable)


** Sinovac
gen sinovac=0 if h16051==3 & h16053==.
replace sinovac=1 if h16051==3 & h16053==3
replace sinovac=2 if h16051==3 & h16053==4
replace sinovac=3 if h16051==3 & h16053==1
replace sinovac=4 if h16051==3 & h16053==6
replace sinovac=5 if h16051==3 & h16053==5
replace sinovac=6 if h16051==3 & h16053==2

label define sinovac 0"not booster" 1"sinovac" 2"aztra" 3"pfizer" 4"moderna" 5"cansino" 6"sputnik", modify
label values sinovac sinovac

** prepare varibale for adjusted model
gen sinovacmodelo=sinovac 
recode sinovacmodelo (0=4) (2=3) (3=5) (4=5) (5=3) (6=3)
***** Recod for colapse heterologous
recode sinovacmodelo (3=5) 
label values sinovacmodelo booster 
tab sinovacmodelo sinovac


** Astra Zeneca
gen aztra=0 if h16051==4 & h16053==.
replace aztra=1 if h16051==4 & h16053==4
replace aztra=2 if h16051==4 & h16053==2
replace aztra=3 if h16051==4 & h16053==1
replace aztra=4 if h16051==4 & h16053==5
replace aztra=5 if h16051==4 & h16053==6
replace aztra=6 if h16051==4 & h16053==3

label define aztra 0"not booster" 1"aztra" 2"sputnik" 3"pfizer" 4"cansino" 5"moderna" 6"sinovac", modify
label values aztra aztra

** prepare varibale for adjusted model
gen aztramodelo=aztra
recode aztramodelo (0=4) (4=2) (5=3) (6=.)
***** Recod for colapse homologous
recode aztramodelo (2=1)
label values aztramodelo booster 
tab aztramodelo aztra


** Cansino
gen cansino=0 if h16051==5 & h16052==.
replace cansino=1 if h16051==5 & h16052==5
replace cansino=2 if h16051==5 & h16052==4
replace cansino=3 if h16051==5 & h16052==6
replace cansino=4 if h16051==5 & h16052==1
replace cansino=5 if h16051==5 & h16052==2
replace cansino=6 if h16051==5 & h16052==3
replace cansino=7 if h16051==5 & h16052==7

label define cansino 0"not booster" 1"cansino" 2"aztra" 3"moderna" 4"pfizer" 5"sputnik" 6"sinovac" 7"JJ", modify
label values cansino cansino

** prepare varibale for adjusted model
gen cansinomodelo=cansino
recode cansinomodelo (0=4) (5=2) (7=2) (4=3) (6=.)
***** Recod for colapse homologous
recode cansinomodelo (2=1)
label values cansinomodelo booster 
tab cansinomodelo cansino

**********************************************
**********************************************

** prepare varibale "Age Group" for models
gen edadcat_nueva=h0303
recode edadcat_nueva 1/17=. 18/59=1 60/100=2
label define edadcat_nueva 1"18 a 59 años" 2"60 años o más"
label values edadcat_nueva edadcat_nueva
tab edadcat_nueva

**********************************************
**********************************************

** Exploration of vaccination dates

**********************************************
**********************************************

** 1st Dose 
** Day
 gen dia1dosis= h16041d
 recode dia1dosis 99=.
** Month
 gen mes1dosis= h16041m
 recode mes1dosis 99=.
** Year
 gen año1dosis=h16041a
 recode año1dosis 9999=.

egen mesaño1 =concat(mes1dosis año1dosis) if mes1dosis!=. & año1dosis!=., punct("-")
tab mesaño1
gen mesaño1_espejo=monthly(mesaño1, "MY")
format mesaño1_espejo %tm 
tab mesaño1_espejo
 
** 2nd 
** Day
 gen dia2dosis= h16042d
 recode dia2dosis 77=. 99=.
** Month
 gen mes2dosis= h16042m
 recode mes2dosis 77=. 99=.
** Year
 gen año2dosis=h16042a
 recode año2dosis 7777=. 9999=.

egen mesaño2 =concat(mes2dosis año2dosis) if mes2dosis!=. & año2dosis!=., punct("-")
tab mesaño2
gen mesaño2_espejo=monthly(mesaño2, "MY")
format mesaño2_espejo %tm 
tab mesaño2_espejo
 
** 3rd Dose
** Day
 gen dia3dosis= h16043d
 recode dia3dosis 77=. 99=.
** Month
 gen mes3dosis= h16043m
 recode mes3dosis 77=. 99=.
** Year
 gen año3dosis=h16043a
 recode año3dosis 7777=. 9999=.

egen mesaño3 =concat(mes3dosis año3dosis) if mes3dosis!=. & año3dosis!=., punct("-")
tab mesaño3
gen mesaño3_espejo=monthly(mesaño3, "MY")
format mesaño3_espejo %tm 
tab mesaño3_espejo


** SAMPLE COLLECTION DATE // Ensanut record 4 visits for sample collection
** Start Date (1st)
tab fecha_ini_1
gen fecha_ini_1_espejo =date(fecha_ini_1, "DMY") 
format fecha_ini_1_espejo %td
tab fecha_ini_1_espejo

** only month and year
gen mesmuestra1=month(fecha_ini_1_espejo) 
tab mesmuestra1
gen añomuestra1=year(fecha_ini_1_espejo) 
tab añomuestra1
tab mesmuestra1 añomuestra1
recode añomuestra1 2021=2022
tab añomuestra1

** 1st sample collection date
egen fechamuestra1 =concat(mesmuestra1 añomuestra) if mesmuestra1!=. & añomuestra1!=., punct("-")
tab fechamuestra1
gen fechamuestra1_espejo=monthly(fechamuestra1, "MY")
format fechamuestra1_espejo %tm 
tab fechamuestra1_espejo

** 2nd date
tab fecha_ini_2 
gen fecha_ini_2_espejo =date(fecha_ini_2, "DMY") 
format fecha_ini_2_espejo %td
tab fecha_ini_2_espejo 

** only month and year
gen mesmuestra2=month(fecha_ini_2_espejo) 
tab mesmuestra2
gen añomuestra2=year(fecha_ini_2_espejo)
tab añomuestra2

** 2nd sample collection date
egen fechamuestra2 =concat(mesmuestra2 añomuestra2) if mesmuestra2!=. & añomuestra2!=., punct("-")
tab fechamuestra2
gen fechamuestra2_espejo=monthly(fechamuestra2, "MY")
format fechamuestra2_espejo %tm 
tab fechamuestra2_espejo

** 3rd date
tab fecha_ini_3 
gen fecha_ini_3_espejo =date(fecha_ini_3, "DMY") 
format fecha_ini_3_espejo %td
tab fecha_ini_3_espejo 

** only month and year
gen mesmuestra3=month(fecha_ini_3_espejo) 
tab mesmuestra3
gen añomuestra3=year(fecha_ini_3_espejo)
tab añomuestra3

** 3rd date sample collection date
egen fechamuestra3 =concat(mesmuestra3 añomuestra3) if mesmuestra3!=. & añomuestra3!=., punct("-")
tab fechamuestra3
 gen fechamuestra3_espejo=monthly(fechamuestra3, "MY")
format fechamuestra3_espejo %tm 
tab fechamuestra3_espejo

** Last date (Fourth)
tab fecha_ini_4
gen fecha_ini_4_espejo =date(fecha_ini_4, "DMY") 
format fecha_ini_4_espejo %td
tab fecha_ini_4_espejo 

** only month and year
gen mesmuestra4=month(fecha_ini_4_espejo) 
tab mesmuestra4
gen añomuestra4=year(fecha_ini_4_espejo)
tab añomuestra4

** Last (4th) date sample collection date
egen fechamuestra4 =concat(mesmuestra4 añomuestra4) if mesmuestra4!=. & añomuestra4!=., punct("-")
tab fechamuestra4
 gen fechamuestra4_espejo=monthly(fechamuestra4, "MY")
format fechamuestra4_espejo %tm 
tab fechamuestra4_espejo


** Final Date // to know the result
tab resultado_1
gen fechafinal1= fechamuestra1_espejo if resultado_1==1
tab fechafinal1
tab fechafinal1,m 
replace fechafinal1 = fechamuestra2_espejo if fechafinal1==. & resultado_2==1

tab resultado_3
replace fechafinal1 = fechamuestra3_espejo if fechafinal1==. & resultado_3==1
tab fechafinal1
tab fechafinal1, m

tab resultado_4
replace fechafinal1 = fechamuestra4_espejo if fechafinal1==. & resultado_4==1
format fechafinal1 %tm 
tab fechafinal1
tab fechafinal1, m

** Difference in months

** fm_f1 -> final date sample - 1st dose date 
gen fm_f1_meses= fechafinal1-mesaño1_espejo
tab fm_f1_meses

** fm_f2 -> final date sample - 2nd dose date
gen fm_f2_meses= fechafinal1-mesaño2_espejo
tab fm_f2_meses

** fm_f3 -> final date sample - 3rd dose date
gen fm_f3_meses= fechafinal1-mesaño3_espejo
tab fm_f3_meses


*************************************
** meses_final // by type of vaccine
*************************************

** CanSino
gen meses_final=.
replace meses_final= fm_f1_meses if cansino==0
** Booster 
replace meses_final= fm_f2_meses if inrange(cansino,1,7)
** Pfzer
** Booster 
replace meses_final= fm_f3_meses if inrange(pfizer,1,7)
** Without booster
replace meses_final= fm_f2_meses if pfizer==0

** Sputnik
** Booster 
replace meses_final= fm_f3_meses if inrange(sputnik,1,7)
** Without Booster
replace meses_final= fm_f2_meses if sputnik==0

** Sinovac
** Booster
replace meses_final= fm_f3_meses if inrange(sinovac,1,7)
** Without Booster
replace meses_final= fm_f2_meses if sinovac==0

** Astra Zeneca
** Booster
replace meses_final= fm_f3_meses if inrange(aztra,1,7)
** Without booster
replace meses_final= fm_f2_meses if aztra==0

tab meses_final
tab meses_final, m


** booster
gen booster=1
replace booster=0 if cansino==0
replace booster=0 if pfizer==0
replace booster=0 if sputnik==0
replace booster=0 if aztra==0
replace booster=0 if sinovac==0
		  
** tipobooster // for imputed
gen tipobooster=h16053
recode tipobooster .=0
replace tipobooster=h16052 if h16051==5
recode tipobooster (.=0)
label define tipobooster 0"not booster" 1"Pfizer" 2"Sputnik" 3"Sinovac" 4"AztraZeneca" 5"CanSino" 6"Moderna" 7"J&J", modify
label values tipobooster tipobooster 
tab tipobooster


**********************************************
**********************************************

** Explore date  

**********************************************
**********************************************

* fecha // without vs with date
gen fecha=0 if meses_final==.
replace fecha=1 if meses_final!=.
tab fecha


**********************************************
**********************************************

** Imputed

**********************************************
**********************************************

ccmatch edadcat entidad municipio tipobooster h16051, cc(fecha)
*Total Matches: 200

sort match
br fecha meses_final match
sort match fecha
br fecha meses_final match
replace meses_final=meses_final[_n+1] if fecha==0 & match==match[_n+1] & match!=. & meses_final==. 
*(200 real changes made) confirm imputed

drop if meses_final==.
*(320 observations deleted)

**********************************************
** Declare database to weight
**********************************************

svyset [pweight =ponde_c_indre], strata(est_sel) psu(upm) singleunit(centered) 

**********************************************
**********************************************


** Table 1. Vaccinated report by sociodemographic characteristics in the Mexican population, ENSANUT 2022.

* (IMPUTED) n=2,953
**********************************************
**********************************************

svy: tab sexo, percent ci
svy: tab edadcat, percent ci
svy: tab escolaridad, percent ci
svy: tab region, percent ci
svy: tab estrato, percent ci
svy: tab nseF, percent ci
svy: tab seguridad, percent ci
svy: tab ocupacion, percent ci



**********************************************
**********************************************
** Table 2.  Geometric Mean Titer anti-S protein antibodies by age group and type of booster in the Mexican population with a primary complete scheme, ENSANUT 2022.
**********************************************
**********************************************


** Age group: Adults and older
*Pfizer
svy: regress var_log i.pfizermodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2) 
margins pfizermodelo, atmeans asbalanced

*Sputnik
svy: regress var_log i.sputnikmodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2)
margins sputnikmodelo, atmeans asbalanced

*Sinovac
svy: regress var_log i.sinovacmodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2)
margins sinovacmodelo, atmeans asbalanced

*Aztra Zeneca
svy: regress var_log i.aztramodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2)
margins aztramodelo, atmeans asbalanced

*CanSino
svy: regress var_log i.cansinomodelo meses_final i.protN h0303  if inrange(edadcat_nueva, 1,2)
margins cansinomodelo, atmeans asbalanced



**********************************************
**********************************************

*** Supplementary material

**********************************************
**********************************************


**********************************************
**********************************************
** Table S1. Frequencies and percentages by sociodemographic characteristics and final date (in months) between primary complete scheme, booster and sample collection date in the Mexican population, ENSANUT 2022. 
**********************************************
**********************************************

tab sexo fecha, col chi2 
tab edadcat fecha, col chi2 
tab escolaridad fecha, col chi2 
tab region fecha, col chi2 
tab estrato fecha, col chi2 
tab nseF fecha, col chi2 
tab seguridad fecha, col chi2 
tab ocupacion fecha, col chi2 


**********************************************
**********************************************
** Table S2. Frequencies and percentages anti-S protein by age group in primary complete scheme, with and without booster of the Mexican population participants in ENSANUT 2022.
**********************************************
**********************************************

tabout pfizer aztra sinovac sputnik cansino edadcat_nueva using tables2.xls, svy percent c(col ci) replace  

tabout pfizer aztra sinovac sputnik cansino edadcat_nueva using tables2.xls, c(freq) append 

**********************************************
** Table S2.1 Frequencies and percentages (models imputed)

tabout pfizermodelo aztramodelo sinovacmodelo sputnikmodelo cansinomodelo edadcat_nueva using tableS2.1model_imputed.xls, svy percent c(col ci) replace  


tabout pfizermodelo aztramodelo sinovacmodelo sputnikmodelo cansinomodelo edadcat_nueva using tableS2.1model_imputed.xls, c(freq) append 

** Table S2.2 Frequencies and percentages 

tabout pfizermodelo aztramodelo sinovacmodelo sputnikmodelo cansinomodelo edadcat_nueva if fecha==1 using tableS2.2modelo.xls, svy percent c(col ci) replace  


tabout pfizermodelo aztramodelo sinovacmodelo sputnikmodelo cansinomodelo edadcat_nueva if fecha==1 using tableS2.2modelo.xls, c(freq) append 


**********************************************
**********************************************
** Table S3. Geometric Mean Titer anti-S protein antibodies by age group and type of booster in the Mexican population with a primary complete scheme, ENSANUT 2022.
**********************************************
**********************************************



** Group: Adults and older
*Pfizer
svy: regress var_log i.pfizermodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2)  & fecha==1
margins pfizermodelo, atmeans asbalanced

*Sputnik
svy: regress var_log i.sputnikmodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2) & fecha==1
margins sputnikmodelo, atmeans asbalanced

*Sinovac
svy: regress var_log i.sinovacmodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2) & fecha==1
margins sinovacmodelo, atmeans asbalanced

*Aztra Zeneca
svy: regress var_log i.aztramodelo meses_final i.protN h0303 if inrange(edadcat_nueva, 1,2) & fecha==1
margins aztramodelo, atmeans asbalanced

*CanSino
svy: regress var_log i.cansinomodelo meses_final i.protN h0303  if inrange(edadcat_nueva, 1,2) & fecha==1
margins cansinomodelo, atmeans asbalanced


**********************************************
**********************************************
** Table S4. Geometric Mean Titer anti-S protein antibodies by age group (Adults=18-59, Older=60 years or more) and type of booster in the Mexican population with a primary complete scheme, ENSANUT 2022.
**********************************************
**********************************************


** Age group 1: 18-59 years
** pfzer
svy: regress var_log i.pfizermodelo meses_final i.protN if edadcat_nueva==1
margins pfizermodelo, atmeans asbalanced

** Sputnik
svy: regress var_log i.sputnikmodelo meses_final i.protN if edadcat_nueva==1
margins sputnikmodelo, atmeans asbalanced

** Sinovac 
svy: regress var_log i.sinovacmodelo meses_final i.protN if edadcat_nueva==1
margins sinovacmodelo, atmeans asbalanced

** Astra Zeneca
svy: regress var_log i.aztramodelo meses_final i.protN if edadcat_nueva==1
margins aztramodelo, atmeans asbalanced

** Cansino 
svy: regress var_log i.cansinomodelo meses_final i.protN if edadcat_nueva==1
margins cansinomodelo, atmeans asbalanced


** Age group 2: 60 years or more

** pfzer
svy: regress var_log i.pfizermodelo meses_final i.protN if edadcat_nueva==2
margins pfizermodelo, atmeans asbalanced

** Sputnik
svy: regress var_log i.sputnikmodelo meses_final i.protN if edadcat_nueva==2
margins sputnikmodelo, atmeans asbalanced

** Sinovac 
svy: regress var_log i.sinovacmodelo meses_final i.protN if edadcat_nueva==2
margins sinovacmodelo, atmeans asbalanced

** Astra Zeneca
svy: regress var_log i.aztramodelo meses_final i.protN if edadcat_nueva==2
margins aztramodelo, atmeans asbalanced

** Cansino 
svy: regress var_log i.cansinomodelo meses_final i.protN if edadcat_nueva==2
margins cansinomodelo, atmeans asbalanced



