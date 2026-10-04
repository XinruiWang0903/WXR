version 17.0
clear all
set more off
set matsize 11000
set seed 20260924

*******************************************************
* 0. 路径与依赖
*******************************************************

* 安装常用扩展命令
capture which winsor2
if _rc ssc install winsor2, replace

capture which ftools
if _rc ssc install ftools, replace

capture which reghdfe
if _rc ssc install reghdfe, replace

capture which esttab
if _rc ssc install estout, replace

capture which ivreg2
if _rc ssc install ivreg2, replace

capture which ranktest
if _rc ssc install ranktest, replace

capture which ivreghdfe
if _rc ssc install ivreghdfe, replace

capture which xsmle
if _rc ssc install xsmle, replace

capture which ddml
if _rc ssc install ddml, replace

capture which pystacked
if _rc ssc install pystacked, replace

*******************************************************
* 1. 数据读取与面板设定
*******************************************************

sort id year
xtset id year
destring CEE GF LES INL DOW GOV FDE TIL GB IS UI CG TI, replace 

label var CEE "能源碳排放效率"
label var GF  "绿色金融"
label var LES "规模经济水平"
label var INL "信息化水平"
label var DOW "对外开放程度"
label var GOV "政府干预力度"
label var FDE "金融发展水平"
label var TIL "科技投入水平"
label var GB  "替代绿色金融指标"
label var IS  "产业绿色转型"
label var UI  "产业结构高级化"
label var CG  "绿色资本引导"
label var TI  "绿色技术创新"
label var Boundary "边界城市"
label var OldIndustrial "老工业基地"

global CTRL "LES INL DOW GOV FDE TIL"

* 样本结构核验
noi di "========================================"
noi di "样本结构检查"
noi di "========================================"
count
noi di "总观测数 = " r(N)

egen tagcity = tag(id)
count if tagcity
noi di "城市数 = " r(N)
drop tagcity

isid id year
tab year
tab Boundary
tab OldIndustrial
tab ResourceType

misstable summarize CEE GF LES INL DOW GOV FDE TIL GB IS UI CG TI

* 主样本校验
count
assert r(N)==2959
quietly levelsof id, local(citylist)
local ncity : word count `citylist'
assert `ncity'==269

save "$OUT/analysis_ready.dta", replace


*******************************************************
* 2. 描述性统计、相关性与多重共线性
*******************************************************

use "$OUT/analysis_ready.dta", clear
xtset id year

* 描述性统计
estpost summarize CEE GF LES INL DOW GOV FDE TIL
esttab using "$OUT/Table3_描述性统计.rtf", replace ///
    cells("count(fmt(0)) mean(fmt(4)) sd(fmt(4)) min(fmt(4)) max(fmt(4))") ///
    nonumber nomtitle noobs label

* 相关系数
pwcorr CEE GF LES INL DOW GOV FDE TIL, sig star(0.05)

estpost correlate CEE GF LES INL DOW GOV FDE TIL, matrix
esttab using "$OUT/附表_相关系数.rtf", replace ///
    unstack not noobs compress

* VIF
reg CEE GF LES INL DOW GOV FDE TIL i.year
estat vif

* 面板信息
xtdescribe


*******************************************************
* 3. Hausman 检验
*******************************************************

use "$OUT/analysis_ready.dta", clear
xtset id year

quietly xtreg CEE GF $CTRL i.year, fe
est store FE

quietly xtreg CEE GF $CTRL i.year, re
est store RE

hausman FE RE, sigmamore


*******************************************************
* 4. 基准回归：城市和年份双固定效应
*******************************************************

eststo clear

eststo b1: reghdfe CEE GF, ///
    absorb(id year) vce(cluster id)

eststo b2: reghdfe CEE GF LES, ///
    absorb(id year) vce(cluster id)

eststo b3: reghdfe CEE GF LES INL, ///
    absorb(id year) vce(cluster id)

eststo b4: reghdfe CEE GF LES INL DOW, ///
    absorb(id year) vce(cluster id)

eststo b5: reghdfe CEE GF LES INL DOW GOV, ///
    absorb(id year) vce(cluster id)

eststo b6: reghdfe CEE GF LES INL DOW GOV FDE, ///
    absorb(id year) vce(cluster id)

eststo b7: reghdfe CEE GF LES INL DOW GOV FDE TIL, ///
    absorb(id year) vce(cluster id)

foreach m in b1 b2 b3 b4 b5 b6 b7 {
    estadd local CityFE "Yes": `m'
    estadd local YearFE "Yes": `m'
}

esttab b1 b2 b3 b4 b5 b6 b7 ///
    using "$OUT/Table4_基准回归.rtf", replace ///
    b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(CityFE YearFE N r2, ///
    labels("City FE" "Year FE" "N" "R-squared")) ///
    label compress nogaps


*******************************************************
* 5. 稳健性检验
*******************************************************

use "$OUT/analysis_ready.dta", clear
xtset id year
eststo clear

* 5.1 替换核心解释变量
eststo r1: reghdfe CEE GB $CTRL, ///
    absorb(id year) vce(cluster id)

* 5.2 剔除直辖市
eststo r2: reghdfe CEE GF $CTRL ///
    if !inlist(city,"北京市","天津市","上海市","重庆市"), ///
    absorb(id year) vce(cluster id)

* 5.3 核心解释变量滞后一期
eststo r3: reghdfe CEE L.GF $CTRL, ///
    absorb(id year) vce(cluster id)

foreach m in r1 r2 r3 {
    estadd local Controls "Yes": `m'
    estadd local CityFE "Yes": `m'
    estadd local YearFE "Yes": `m'
}


*******************************************************
* 6. 内生性检验：工具变量法
*******************************************************

* 以绿色金融一期滞后项为工具变量
eststo r4: ivreghdfe CEE $CTRL (GF=L.GF), ///
    absorb(id year) cluster(id) first

estadd local Controls "Yes": r4
estadd local CityFE "Yes": r4
estadd local YearFE "Yes": r4

esttab r1 r2 r3 r4 ///
    using "$OUT/Table5_稳健性与内生性.rtf", replace ///
    b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(Controls CityFE YearFE N, ///
    labels("Controls" "City FE" "Year FE" "N")) ///
    label compress nogaps

* 单独输出第一阶段及识别统计量
ivreghdfe CEE $CTRL (GF=L.GF), ///
    absorb(id year) cluster(id) first
	

*******************************************************
* 7. 双重机器学习 DML
*******************************************************

capture noisily {

    use "$OUT/analysis_ready.dta", clear
    xtset id year

    * 先剔除城市和年份固定效应
    foreach v in CEE GF LES INL DOW GOV FDE TIL {
        quietly reghdfe `v', absorb(id year)
        predict double R_`v', resid
    }

    * 控制变量二次项
    foreach v in LES INL DOW GOV FDE TIL {
        gen double R_`v'_sq = R_`v'^2
    }

    global X1 "R_LES R_INL R_DOW R_GOV R_FDE R_TIL"
    global X2 "R_LES R_INL R_DOW R_GOV R_FDE R_TIL R_LES_sq R_INL_sq R_DOW_sq R_GOV_sq R_FDE_sq R_TIL_sq"

    tempname dmlmem
    postfile `dmlmem' str30 model double coef se p ///
        using "$OUT/Table6_DML_results.dta", replace

    * 7.1 随机森林：5折，一次项
    ddml init partial, kfolds(5)
    ddml E[Y|X]: pystacked R_CEE $X1, type(reg) methods(rf)
    ddml E[D|X]: pystacked R_GF  $X1, type(reg) methods(rf)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("RF_5fold_linear") (bb) (ss) (pp)

    * 7.2 随机森林：5折，二次项
    ddml clear
    ddml init partial, kfolds(5)
    ddml E[Y|X]: pystacked R_CEE $X2, type(reg) methods(rf)
    ddml E[D|X]: pystacked R_GF  $X2, type(reg) methods(rf)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("RF_5fold_quadratic") (bb) (ss) (pp)

    * 7.3 随机森林：8折，一次项
    ddml clear
    ddml init partial, kfolds(8)
    ddml E[Y|X]: pystacked R_CEE $X1, type(reg) methods(rf)
    ddml E[D|X]: pystacked R_GF  $X1, type(reg) methods(rf)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("RF_8fold_linear") (bb) (ss) (pp)

    * 7.4 随机森林：8折，二次项
    ddml clear
    ddml init partial, kfolds(8)
    ddml E[Y|X]: pystacked R_CEE $X2, type(reg) methods(rf)
    ddml E[D|X]: pystacked R_GF  $X2, type(reg) methods(rf)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("RF_8fold_quadratic") (bb) (ss) (pp)

    * 7.5 ElasticNet
    ddml clear
    ddml init partial, kfolds(5)
    ddml E[Y|X]: pystacked R_CEE $X2, type(reg) methods(elasticcv)
    ddml E[D|X]: pystacked R_GF  $X2, type(reg) methods(elasticcv)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("ElasticNet") (bb) (ss) (pp)

    * 7.6 Gradient Boosting
    ddml clear
    ddml init partial, kfolds(5)
    ddml E[Y|X]: pystacked R_CEE $X2, type(reg) methods(gradboost)
    ddml E[D|X]: pystacked R_GF  $X2, type(reg) methods(gradboost)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("GradientBoost") (bb) (ss) (pp)

    * 7.7 LassoCV
    ddml clear
    ddml init partial, kfolds(5)
    ddml E[Y|X]: pystacked R_CEE $X2, type(reg) methods(lassocv)
    ddml E[D|X]: pystacked R_GF  $X2, type(reg) methods(lassocv)
    ddml crossfit
    quietly ddml estimate, robust
    matrix B=e(b)
    matrix V=e(V)
    scalar bb=B[1,1]
    scalar ss=sqrt(V[1,1])
    scalar pp=2*normal(-abs(bb/ss))
    post `dmlmem' ("LassoCV") (bb) (ss) (pp)

    postclose `dmlmem'

    use "$OUT/Table6_DML_results.dta", clear
    gen sig=""
    replace sig="***" if p<0.01
    replace sig="**"  if p>=0.01 & p<0.05
    replace sig="*"   if p>=0.05 & p<0.10

    export excel using "$OUT/Table6_DML_results.xlsx", ///
        firstrow(variables) replace

    export delimited using "$OUT/Table6_DML_results.csv", replace
}


*******************************************************
* 8. 机制分析
*******************************************************

use "$OUT/analysis_ready.dta", clear
xtset id year
eststo clear

* 8.1 产业绿色转型
eststo m1: reghdfe IS GF $CTRL, ///
    absorb(id year) vce(cluster id)

eststo m2: reghdfe CEE GF IS $CTRL, ///
    absorb(id year) vce(cluster id)

* 8.2 产业结构高级化
eststo m3: reghdfe UI GF $CTRL, ///
    absorb(id year) vce(cluster id)

eststo m4: reghdfe CEE GF UI $CTRL, ///
    absorb(id year) vce(cluster id)

* 8.3 绿色资本引导
eststo m5: reghdfe CG GF $CTRL, ///
    absorb(id year) vce(cluster id)

eststo m6: reghdfe CEE GF CG $CTRL, ///
    absorb(id year) vce(cluster id)

* 8.4 绿色技术创新
eststo m7: reghdfe TI GF $CTRL, ///
    absorb(id year) vce(cluster id)

eststo m8: reghdfe CEE GF TI $CTRL, ///
    absorb(id year) vce(cluster id)

foreach m in m1 m2 m3 m4 m5 m6 m7 m8 {
    estadd local Controls "Yes": `m'
    estadd local CityFE "Yes": `m'
    estadd local YearFE "Yes": `m'
}

esttab m1 m2 m3 m4 m5 m6 m7 m8 ///
    using "$OUT/Table7_机制分析.rtf", replace ///
    b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(Controls CityFE YearFE N r2, ///
    labels("Controls" "City FE" "Year FE" "N" "R-squared")) ///
    label compress nogaps


*******************************************************
* 9. 异质性分析
*******************************************************

use "$OUT/analysis_ready.dta", clear
xtset id year
eststo clear

* 9.1 边界城市 / 非边界城市
eststo h1: reghdfe CEE GF $CTRL ///
    if Boundary==1, absorb(id year) vce(cluster id)

eststo h2: reghdfe CEE GF $CTRL ///
    if Boundary==0, absorb(id year) vce(cluster id)

* 9.2 老工业基地 / 非老工业基地
eststo h3: reghdfe CEE GF $CTRL ///
    if OldIndustrial==1, absorb(id year) vce(cluster id)

eststo h4: reghdfe CEE GF $CTRL ///
    if OldIndustrial==0, absorb(id year) vce(cluster id)

* 9.3 四类资源型城市
eststo h5: reghdfe CEE GF $CTRL ///
    if ResourceType=="成熟型", absorb(id year) vce(cluster id)

eststo h6: reghdfe CEE GF $CTRL ///
    if ResourceType=="成长型", absorb(id year) vce(cluster id)

eststo h7: reghdfe CEE GF $CTRL ///
    if ResourceType=="再生型", absorb(id year) vce(cluster id)

eststo h8: reghdfe CEE GF $CTRL ///
    if ResourceType=="衰退型", absorb(id year) vce(cluster id)

foreach m in h1 h2 h3 h4 h5 h6 h7 h8 {
    estadd local Controls "Yes": `m'
    estadd local CityFE "Yes": `m'
    estadd local YearFE "Yes": `m'
}

esttab h1 h2 h3 h4 h5 h6 h7 h8 ///
    using "$OUT/Table8_异质性分析.rtf", replace ///
    b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    stats(Controls CityFE YearFE N r2, ///
    labels("Controls" "City FE" "Year FE" "N" "R-squared")) ///
    label compress nogaps


*******************************************************
* 10. 空间矩阵导入
*******************************************************

* 邻接矩阵
clear
import delimited using "$DATA/空间邻接权重矩阵.csv", ///
    varnames(1) encoding(utf8)
drop id
mkmat c*, matrix(W_adj)

* 反距离矩阵
clear
import delimited using "$DATA/反距离地理矩阵.csv", ///
    varnames(1) encoding(utf8)
drop id
mkmat c*, matrix(W_inv)


*******************************************************
* 11. Global Moran's I
*******************************************************

use "$OUT/analysis_ready.dta", clear
sort id year
xtset id year

tempname moranmem
postfile `moranmem' int year double Moran_adj Moran_inv ///
    using "$OUT/Moran_I_by_year.dta", replace

forvalues yy=2013/2023 {

    preserve

    keep if year==`yy'
    sort id
    mkmat CEE, matrix(Y)

    mata:
        y  = st_matrix("Y")
        Wa = st_matrix("W_adj")
        Wi = st_matrix("W_inv")
        e  = y :- mean(y)
        n  = rows(y)

        Ia = (n/sum(Wa))*((e'*Wa*e)/(e'*e))
        Ii = (n/sum(Wi))*((e'*Wi*e)/(e'*e))

        st_numscalar("MIa",Ia)
        st_numscalar("MIi",Ii)
    end

    post `moranmem' (`yy') (MIa) (MIi)

    restore
}

postclose `moranmem'

preserve
use "$OUT/Moran_I_by_year.dta", clear
list, clean noobs
export excel using "$OUT/Moran_I_by_year.xlsx", ///
    firstrow(variables) replace
restore


*******************************************************
* 12. 空间模型选择与空间杜宾模型
*******************************************************

use "$OUT/analysis_ready.dta", clear
sort id year
xtset id year

* 普通FE / RE Hausman
quietly xtreg CEE GF $CTRL i.year, fe
est store FE_sp

quietly xtreg CEE GF $CTRL i.year, re
est store RE_sp

hausman FE_sp RE_sp, sigmamore

* 邻接矩阵：SAR / SEM / SDM
capture noisily xsmle CEE GF $CTRL, ///
    wmat(W_adj) model(sar) fe type(both)
capture estimates store SAR_adj

capture noisily xsmle CEE GF $CTRL, ///
    wmat(W_adj) model(sem) fe type(both)
capture estimates store SEM_adj

capture noisily xsmle CEE GF $CTRL, ///
    wmat(W_adj) model(sdm) fe type(both)
capture estimates store SDM_adj0

capture noisily lrtest SDM_adj0 SAR_adj
capture noisily lrtest SDM_adj0 SEM_adj


*******************************************************
* 13. 空间杜宾模型：两种权重矩阵
*******************************************************

* 先做空间滞后方向核验
eststo clear

eststo sx1: reghdfe CEE GF WGF_adj $CTRL, ///
    absorb(id year) vce(cluster id)

eststo sx2: reghdfe CEE GF WGF_inv $CTRL, ///
    absorb(id year) vce(cluster id)

esttab sx1 sx2 ///
    using "$OUT/Table9A_空间滞后方向核验.rtf", replace ///
    b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    label compress

* 13.1 空间邻接矩阵
xsmle CEE GF $CTRL, ///
    wmat(W_adj) model(sdm) fe type(both) effects

estimates store SDM_adj
estimates save "$OUT/SDM_adj", replace

* 13.2 反距离地理矩阵
xsmle CEE GF $CTRL, ///
    wmat(W_inv) model(sdm) fe type(both) effects

estimates store SDM_inv
estimates save "$OUT/SDM_inv", replace

* 输出主方程
esttab SDM_adj SDM_inv ///
    using "$OUT/Table9B_SDM主方程.rtf", replace ///
    b(%9.3f) se(%9.3f) ///
    star(* 0.10 ** 0.05 *** 0.01) ///
    label compress

	
*******************************************************
* 14. 100—1500 km 距离阈值分析
*******************************************************

use "$OUT/analysis_ready.dta", clear
xtset id year

tempname distmem
postfile `distmem' int distance double coef se p ///
    using "$OUT/distance_threshold_results.dta", replace

forvalues d=100(100)1500 {

    quietly reghdfe CEE GF WGF_`d' $CTRL, ///
        absorb(id year) vce(cluster id)

    scalar bb = _b[WGF_`d']
    scalar ss = _se[WGF_`d']
    scalar pp = 2*ttail(e(df_r),abs(bb/ss))

    post `distmem' (`d') (bb) (ss) (pp)
}

postclose `distmem'

use "$OUT/distance_threshold_results.dta", clear

gen upper = coef + 1.96*se
gen lower = coef - 1.96*se

gen sig = ""
replace sig="***" if p<0.01
replace sig="**"  if p>=0.01 & p<0.05
replace sig="*"   if p>=0.05 & p<0.10

list, clean noobs

export excel using "$OUT/Figure2_距离阈值数据.xlsx", ///
    firstrow(variables) replace

export delimited using "$OUT/Figure2_距离阈值数据.csv", replace

twoway ///
    (rarea upper lower distance, color(gs13%55) lcolor(none)) ///
    (connected coef distance, msymbol(circle) lwidth(medthick)), ///
    yline(0, lcolor(black) lpattern(dash)) ///
    xtitle("Geographical distance threshold (km)") ///
    ytitle("Spatial spillover coefficient") ///
    xlabel(100(200)1500) ///
    legend(order(1 "95% confidence interval" 2 "Spatial spillover effect")) ///
    graphregion(color(white))

graph export "$OUT/Figure2_距离阈值.png", ///
    replace width(2400)

graph save "$OUT/Figure2_距离阈值.gph", replace