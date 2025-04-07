import time
import random as rd
seed =int(rd.random()*300)
with open("mRNApoly_ninfo_1bp.inp") as ff:
    ff_=ff.readlines()
    ff_[589]="n_seed = {}\n".format(seed)
with open("mRNApoly_ninfo_1bp.inp", mode='w') as gg:
    gg.writelines(ff_)
### 上記でrandom seed を変えた。    
import subprocess
import time
import tqdm
for i in tqdm.tqdm(range(1,120)):
    subprocess.run("/home/yamada/cafemol/bin/cafemol mRNApoly_ninfo_{}bp.inp".format(i), shell=True)
    time.sleep(20)
    with open("mRNApoly_ninfo_{}bp.inp".format(i)) as f:##i+1番目のinpファイルを作るためi番目のinpを参照
        l = f.readlines()###lはあるinpファイルに関してリスト化したもので、inpファイルそのものとは別
        l[5]="filename = mr_output_ninf_{}\n".format(i+1)###outputファイル名の編集
        l[7]="path_natinfo = ./RNApoly_ninfo_edit/一時逆方向/{}bp_move\n".format(i+1)####ninfoファイルを変える 逆方向(正しい###方向)か順方向(誤った方向)によって変える。
        l[546]="1-16 mr_output_ninf_{}.pdb\n".format(i)###initialstructを新しくできたタンパクを用いる
    time.sleep(10)    
    with open("mRNApoly_ninfo_{}bp.inp".format(i+1), mode='w') as g:
        g.writelines(l)
    time.sleep(10)    
    with open ("mr_output_ninf_{}.pdb".format(i)) as h:###mode="w"は絶対ダメ
        h_ = h.readlines()
        h__ = h_[4753:]###ここからh__の最後五行を編集する
        time.sleep(10)
        for p in range(4746,4751):###4746行目から4750行目まで　原子番号4715から4719まで h__[p]の書き換え なぜかこの部分でバ###グが起こるから
            xx=h__[4743].split()
            x=h__[4744].split()##原子番号4713
            y=h__[4745].split()##原子番号4714
            xx1,xx2,xx3 = float(xx[6]),float(xx[7]),float(xx[8])
            x1,x2,x3 = float(x[6]),float(x[7]),float(x[8])
            y1,y2,y3 = float(y[6]),float(y[7]),float(y[8])
            p_=p-4746###0,1,2,3,4
            z=h__[p].split()
            #z[6],z[7],z[8]=str((x1+y1)/2+p_*2),str((x2+y2)/2+p_*2),str((x3+y3)/2+p_*2)
            w=''
            for q in range(len(z)):
                if q==0:
                    w+=z[q]
                elif q==1:
                    w+=z[q].rjust(7)
                elif q==2:
                    w+=z[q].rjust(4)
                elif q==3:
                    w+=z[q].rjust(5)
                elif q==4:
                    w+=z[q].rjust(2)
                elif q==5:
                    w+=z[q].rjust(4)
                elif q==6:
                    w+=str(round((y1+(y1-x1)*(p_+1)),2)).rjust(12)###小数点第二位まで
                elif q==7:
                    w+=str(round((y2+(y2-x2)*(p_+1)),2)).rjust(8)
                elif q==8:
                    w+=str(round((y3+(y3-x3)*(p_+1)),2)).rjust(8)    
            h__[p]=w+"\n"
            #h__[p]　= h__[p-5]
    time.sleep(10)        
    with open ("mr_output_ninf_{}.pdb".format(i), mode='w') as j:
        j.writelines(h__)
    time.sleep(10)    
    
