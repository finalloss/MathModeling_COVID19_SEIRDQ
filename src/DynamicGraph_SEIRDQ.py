#-*- coding: utf-8 -*-
#导入库
from pylab import *
import math
import numpy as np
import matplotlib.pyplot as plt
#基本参数
Incubation = 7           #感染到住院天数
observing_time = 14      #住院观察天数
Quaranitine_time = 14    #隔离天数
ri = 0.3              #感染范围
I_Infectious_Rate = 0.7#感染概率   
Move_Rate = 0.3       #人员流动参数

Mild_Rate = 0.96       #轻症比例
Severe_Rate = 1 - Mild_Rate   #重症比例
Mild_Lethallty_Rate = 0.0001     #轻症每日死亡率
Severe_Lethallty_Rate = 0.01    #重症每日死亡率

persons = 2000       #社区总人口
days = 30         #历史最大天数
Bed_Rate = 0.15
bed_num = persons*Bed_Rate        #医院总床位
Vaccine_Coverage = 0            #疫苗覆盖率
Ifhalf = 0           #如果医院以某一比例以上床位被占满，发布封城禁令，不再move
Find_Rate = 1/4      #床位被占满比例
IfPierced = 0        #医院是否被击穿，医院床位被占满后不再接纳新患者，直到存在空床位
Iffirst = 1


Superinfection = 0.9        #抗体有效率
Accuracy_Nucleic_acid_testing = 0.5  #核酸检测有效率
Hospital_Recovery = 0.2    #医院治疗能力（单日出院概率）
self_Recovery = 0.05       #个人自治能力

#容器
community = []        #社区人口
latent = []           #潜伏期人群
hospital = []         #住院人群
death = []            #死亡人群
people = []           #总人口

#统计数据
infected_person_byday = []  #每日累计感染人数
admitted_person_byday = []  #每日现存在院人数
dead_person_byday = []      #每日累计死亡人数
cured_person_byday = []     #每日累计治愈人数
Quaranitine_person_byday = []#每日累计隔离人数
infected_person_today = []  #当日感染人数
dead_person_today = []      #当日死亡人数
Quaranitine_person_today = []#每日隔离人数

#position类作为保存person坐标的辅助类
class position:
    x=0
    y=0
    def __init__(self):
        pass
    def setposition(self,x0,y0):
        self.x=x0
        self.y=y0


#person类
class person:
    p = position() #位置坐标
    s = 1     #是否易感，1=易感，0=免疫
    i = 0     #是否患病/传染，1=患病，0=未患病
    c = 0     #是否重症，1=重症，0=轻症，只有在i=1时有效
    d = 0     #是否死亡，1=死亡，0=存活，失去传染能力以及易感性
    Infected_day = -1   #感染日
    Hospital_day = -1   #住院日
    Hospital_time = -1  #留院时长
    Quaranitine_day = -1#隔离日期
    Quaranitine_Mild_Rate = 0.9    #封城后由于物质缺乏轻症比例下降，再进行一次重症判断
    index = 0  #索引号，唯一确定一个人
    Severe_Lethallty_Rate = 0.01
    I_Infectious_Rate = 0.9 
    ri = 0.2
    Quaranitine_first = 0
    #构造方法，idx为索引号，x0、y0为位置横纵坐标
    def __init__(self,idx,x0,y0):
        self.p = position()
        self.index = idx
        self.p.setposition(x0,y0)
    #获得横坐标
    def getx(self):
        return self.p.x
    #获得纵坐标
    def gety(self):
        return self.p.y
    #人员流动
    def move(self):
        self.p.x = self.p.x + Move_Rate * ((1-(-1)) * np.random.random() + (-1))
        self.p.y = self.p.y + Move_Rate * ((1-(-1)) * np.random.random() + (-1))
    #注射疫苗
    def GetVaccined(self):
        self.s = 0
        if np.random.rand() > Vaccine_Coverage:#有概率接种疫苗，代表疫苗接种率
            self.s = 1
        if np.random.rand() > Superinfection:#疫苗有概率失效
            self.s = 1   
    #核酸检测    
    def Nucleic_acid_testing(self,day):
        if self.i == 1 and np.random.rand() > Accuracy_Nucleic_acid_testing:#Accuracy_Nucleic_acid_testing概率检测出阳性
            return 1        #阳性被检测出
        else:
            return 0        #检测为阴性
    #求感染，未感染者计算与潜伏期患者的距离，若小于最大感染距离则有概率被感染
    def If_infected(self,day):
        if self.s == 1:
            for item in latent:
                r = caldistance(self,item)
                if r < self.ri:
                    self.infected(day)                    
    #被感染
    def infected(self,day):
        if self.i == 0:
            if np.random.rand() < self.I_Infectious_Rate:
                self.i = 1
                self.s = 0
                self.Infected_day = day
                self.severity(day)
                infected_person_byday[day] += 1      #统计变量
                latent.append(self)

    #病情分级：轻症、重症
    def severity(self,day):
        IfSeverity = np.random.rand()
        if IfSeverity < Mild_Rate:
            self.c = 0
        else:
            self.c = 1
        self.If_die(day)
    #是否死亡
    def If_die(self,day):
        if self.c==0:  #c为0是轻症
            if np.random.rand() < Mild_Lethallty_Rate:
                self.die(day)
        else:
            if np.random.rand() < self.Severe_Lethallty_Rate:
                self.die(day)
    #住院
    def hospital(self,day):
        self.Hospital_day = day
        admitted_person_byday[day] += 1      #统计变量
        hospital.append(self)
        latent.remove(self)
        community.remove(self)

    #治愈出院，获得免疫
    def cured(self,day):
        self.i = 0
        if np.random.rand() > Superinfection:#若抗体消失，则会再变为易感人群
            self.s = 1
        cured_person_byday[day] += 1         #统计变量        
        if self in hospital:
            self.Hospital_day = -1
            self.Hospital_time = -1
            community.append(self)
            hospital.remove(self)
            admitted_person_byday[day] -= 1      #统计变量
    #死亡
    def die(self,day):
        self.d = 1
        self.s = 0
        self.i = 0
        dead_person_byday[day] += 1          #统计变量        
        death.append(self)
        if self in hospital:
            hospital.remove(self)
            admitted_person_byday[day] -= 1      #统计变量
    #生命周期函数
    def Lifecycle(self,day):
        if(day == 1):
            self.GetVaccined()
        if(day == 10):
            self.I_Infectious_Rate = self.I_Infectious_Rate + 0.2
            self.ri = 2 * self.ri
            self.Severe_Lethallty_Rate = 2 * self.Severe_Lethallty_Rate
        if self.i == 1 and Ifhalf == 0:   #被感染,核酸检测未开始
            d = day-self.Infected_day             
            if d < Incubation:   #潜伏期,人自由移动
                self.move()
            elif d == Incubation:   #住院日，潜伏期过
                self.hospital(day)
            elif d < observing_time+Incubation:    #在住院
                if np.random.random() < Hospital_Recovery:  #单日是否被治愈
                    self.cured(day)
                self.If_die(day)
            elif d == observing_time + Incubation:    #出院日，被治愈
                self.cured(day)
            else:
                pass
        elif self.i == 1 and Ifhalf == 1 and IfPierced == 0: #被感染,核酸检测开始,医疗系统未被击穿，也可能恢复不封城状态
            self.I_Infectious_Rate = self.I_Infectious_Rate - 0.2    #由于发现疫情，人们选择戴口罩，降低传染概率
            if self.Quaranitine_day != -1:  #在解封后的隔离点，移出隔离点
                self.Quaranitine_day == -1
                self.p.x = 0 + ((2-(-2)) * np.random.random() + (-2))
                self.p.y = 0 + ((2-(-2)) * np.random.random() + (-2))   
            IfPositive = self.Nucleic_acid_testing(day)     #进行第day天的核酸检测
            if IfPositive == 0:
                self.move()
            else:    #阳性，住院
                if self not in hospital: #判断是否已经入住
                    self.hospital(day)#未入住医院则进入医院
                    return
                self.Hospital_time = day-self.Hospital_day
                if self.Hospital_time < observing_time:                    
                    self.Severe_Lethallty_Rate = 0.001
                    self.If_die(day)    #住院中
                    if np.random.random() < Hospital_Recovery: #单日是否被治愈
                        self.cured(day)
                elif self.Hospital_time == observing_time:
                    self.cured(day)         #治愈（观察期结束），到时间也会治愈

        elif self.i == 1 and IfPierced == 1:  #被感染,核酸检测开始,医疗系统被击穿，封城不允许出行
            if self not in hospital:       #封城时不在医院
                IfPositive = self.Nucleic_acid_testing(day)#进行第day天的核酸检
                if IfPositive == 0 and self.Quaranitine_day == -1:#未隔离
                    self.Severe_Lethallty_Rate = 0.02 #居家
                    if self.c == 0:#轻症
                        IfSeverity = np.random.rand()
                        if IfSeverity > self.Quaranitine_Mild_Rate:
                            self.c = 1
                    self.If_die(day) #判断是否死亡
                else:#检测出阳性
                    #进入隔离点
                    self.p.x = -15 + ((1-(-1)) * np.random.random() + (-1))
                    self.p.y = -15 + ((1-(-1)) * np.random.random() + (-1))
                    if self.Quaranitine_day == -1:
                        self.Quaranitine_day = day
                    self.Severe_Lethallty_Rate = 0.04
                    if self.Quaranitine_first == 0:
                        Quaranitine_person_byday[day] += 1
                        self.Quaranitine_first = 1
                    if self.c == 0:#轻症
                        IfSeverity = np.random.rand()
                        if IfSeverity > self.Quaranitine_Mild_Rate:
                            self.c = 1
                    self.If_die(day)            #判断是否死亡
                    if np.random.random() < self_Recovery:  #单日是否自愈
                        self.cured(day)
                        self.Quaranitine_day == -1 #离开隔离点
                        self.p.x = ((2-(-2)) * np.random.random() + (-2)) + ((2-(-2)) * np.random.random() + (-2))
                        self.p.y = ((2-(-2)) * np.random.random() + (-2)) + ((2-(-2)) * np.random.random() + (-2))   
            else: #封城时在医院
                self.Hospital_time = day - self.Hospital_day
                if self.Hospital_time < observing_time:
                    self.Severe_Lethallty_Rate = 0.001
                    self.If_die(day)    #住院中
                    if np.random.random() < Hospital_Recovery:  #单日是否被治愈
                        self.cured(day)
                elif self.Hospital_time == observing_time:
                    self.cured(day)         #治愈（观察期结束）
        else:   #若没有被感染
            if IfPierced == 0 and self.d == 0: #未封城
                self.move()
            self.If_infected(day)

#计算两人之间的距离，参数为两person类实例，返回两者距离
def caldistance(p1,p2):
    x1=p1.p.x
    y1=p1.p.y
    x2=p2.p.x
    y2=p2.p.y
    dis = math.sqrt((x1-x2)**2+(y1-y2)**2)
    return dis


def plots():
    X=[]        #横坐标列表
    Y=[]        #纵坐标列表
    C=[]        #散点颜色列表
    Xp = np.arange(0,days)    #统计图横坐标
    for i in range(len(people)):
        X.append(people[i].getx())
        Y.append(people[i].gety())
        if people[i].i==0:
            if people[i].s==1:
                C.append('b')   #未感染者：蓝色 i=0 s=1
            else:
                if people[i].d == 1:
                    C.append('k')   #死者：黑色 i=0 s=0 d=1
                else:
                    C.append('y')   #治愈者：黄色 i=0 s=0 d=0
        else:
            C.append('r')       #潜伏期：红色 i=1
    plt.ion()       #打开动态模式
    plt.clf()       #清除图像
    #支持中文   
    mpl.rcParams['font.sans-serif'] = ['SimHei']
    #用来正常显示负号
    plt.rcParams['axes.unicode_minus'] = False 
    #绘制散点图
    apen=plt.subplot(2,4,1)
    apen.set_title("Day %d" %(day))
    apen.scatter(X,Y,c=C,s=1)
    #绘制累计感染人数图
    bpen=plt.subplot(2,4,2)
    bpen.set_title("累计感染 %.2f%%" %(infected_person_byday[day]/persons*100))
    bpen.plot(Xp[0:day+1],infected_person_byday)
    #绘制当日感染人数图
    bpen=plt.subplot(2,4,3)
    bpen.set_title("每日新增人数 %.2f" %(infected_person_byday[day]-infected_person_byday[day-1]))
    bpen.plot(infected_person_today)
    #绘制实时在院人数图
    cpen=plt.subplot(2,4,4)
    cpen.set_title("实时在院人数 %.2f%%" %(admitted_person_byday[day]/persons*100))
    cpen.plot(Xp[0:day+1],admitted_person_byday)
    #绘制累计死亡人数图
    dpen=plt.subplot(2,4,5)
    dpen.set_title("累计死亡人数 %.2f%%" %(dead_person_byday[day]/persons*100))
    dpen.plot(Xp[0:day+1],dead_person_byday)
    #绘制每日死亡人数图
    epen=plt.subplot(2,4,6)
    epen.set_title("每日死亡 %.2f" %(dead_person_byday[day]-dead_person_byday[day-1]))
    epen.plot(dead_person_today)
    #绘制累计死亡人数图
    fpen=plt.subplot(2,4,7)
    fpen.set_title("累计新增隔离点收纳人数 %.2f%%" %(Quaranitine_person_byday[day]/persons*100))
    fpen.plot(Xp[0:day+1],Quaranitine_person_byday)
    #绘制每日隔离人数图
    gpen=plt.subplot(2,4,8)
    gpen.set_title("每天新增隔离点收纳人数 %.2f" %(Quaranitine_person_byday[day] - Quaranitine_person_byday[day-1]))
    gpen.plot(Quaranitine_person_today)
    plt.tight_layout() 
    plt.pause(0.1)    #暂停
    plt.show()      #显示图像
    #plt.savefig('Day %d.png'%(day))     #保存图像

#产生社区，横纵坐标分布为标准正态分布扩大3倍
Xc=np.random.randn(int(persons/2))*2 
Yc=np.random.randn(int(persons/2))*2
Xd=np.random.randn(int(persons/2))*2 
Yd=np.random.randn(int(persons/2))*2
for i in range(int(persons/2)):
    temp = person(i,Xc[i]+3,Yc[i]+3)
    community.append(temp)
    people.append(temp)
    temp = person(i,Xd[i]-3,Yd[i]-3)
    community.append(temp)
    people.append(temp)
#时间开始
day = 0
#初始化统计变量
cured_person_byday.append(0)
infected_person_byday.append(0)
admitted_person_byday.append(0)
dead_person_byday.append(0)
Quaranitine_person_byday.append(0)
#随机产生零号感染者
community[int(np.random.rand()*persons)].infected(day)
#community[int(np.random.rand()*persons)].infected(day)
print("day0:病毒已选择好零号病人")
#绘制Day 0图像
plots()
#流行开始
#la=1    #潜伏期人数
while day <= days:
    day += 1
    #统计变量列表增加新增1天
    cured_person_byday.append(cured_person_byday[day-1])
    infected_person_byday.append(infected_person_byday[day-1])
    admitted_person_byday.append(admitted_person_byday[day-1])
    dead_person_byday.append(dead_person_byday[day-1])
    Quaranitine_person_byday.append(Quaranitine_person_byday[day-1])
    if admitted_person_byday[day-1] > bed_num*Find_Rate and Iffirst == 1:
        print("day:%d 发现一轻微疾病正在传播，每日核酸检测开始" %day)
        Ifhalf = 1
        Iffirst = 0
    if admitted_person_byday[day-1] > bed_num:
        print("day:%d 警告，医疗系统已被击穿，开始封城" %day)
        IfPierced = 1
    if admitted_person_byday[day-1] < bed_num*Find_Rate and Iffirst == 0:
        print("day:%d 疫情得到控制，free!!!"%day)
        IfPierced = 0
 
    #先处理在院人群
    for i in range(-len(hospital)+1,1):    
        hospital[-i].Lifecycle(day)
    #再处理社区人群
    for i in range(-len(community)+1,1):   
        community[-i].Lifecycle(day)
    infected_person_today.append(infected_person_byday[day] - infected_person_byday[day-1])
    dead_person_today.append(dead_person_byday[day] - dead_person_byday[day-1])
    Quaranitine_person_today.append(Quaranitine_person_byday[day] - Quaranitine_person_byday[day-1])
    #print((infected_person_byday[day]-infected_person_byday[day-1])/la)
    #la=len(latent)
    #画图
    if day == days:
        break
    plots()
print("END")