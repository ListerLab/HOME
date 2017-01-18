import pandas as pd 

import numpy as np
import math
from scipy import stats
from sklearn import preprocessing
from sklearn.externals import joblib
import statsmodels.stats.proportion as sm
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def format_allc(df,classes):
   if classes=="CG":
        filter_col = [col for col in list(df) if col.startswith(('chr',"strand",'pos','mc','h'))]
        df=df[filter_col]
        v=df.chr[0] 
        df1=df[df.strand=='-']
        df2=df[df.strand=='+']
        df1_modpos=df1.pos-1
        df1.loc[:,'pos']=df1_modpos
        df2=df2.drop(["chr","strand"],1)
        df1=df1.drop(["chr","strand"],1)
        df_mod=pd.concat([df2,df1])
        df_mod=df_mod.sort_values(['pos'])
        df=df_mod.groupby(by=['pos']).sum().reset_index()
        df.insert(0,'chr',v)
        filter_col = [col for col in list(df) if col.startswith(('h'))]
        df=df[(df[filter_col]> 2).all(axis=1)]
        df=df.reset_index(drop=True)
   elif classes=="CHN" or classes=="CHG" or classes=="CHH" or classes=="CNN":
        filter_col = [col for col in list(df) if col.startswith(('h'))]
        df=df[(df[filter_col]> 2).all(axis=1)]
        df=df.reset_index(drop=True)
        filter_col = [col for col in list(df) if col.startswith(('chr','pos','mc','h'))]
        df=df[filter_col]
        df=df.reset_index(drop=True)        
      
   return df

def pval_cal_withoutrep(df):
    df=df.rename(columns = {'mc_case_rep1':'mc_case','mc_cont_rep1':'mc_cont','h_case_rep1':'h_case','h_cont_rep1':'h_cont'})
    meth_case=df.mc_case.divide(df.h_case)
    meth_cont=df.mc_cont.divide(df.h_cont)
    
    df['meth_diff']=meth_case.subtract(meth_cont)
    df['meth_case']=meth_case
    df['meth_cont']=meth_cont
    
    df['ztestpval'] =df.apply(lambda r: sm.proportions_ztest(np.array([r.mc_cont,r.mc_case]), np.array([r.h_cont,r.h_case]), value=0, alternative='two-sided')[1], axis=1)
    #df['Fisherpval'] = df.apply(lambda r: stats.fisher_exact([[r.mc, (r.h-r.mc)],[r.mc1,(r.h1-r.mc1)]])[1], axis=1)
    df=df.fillna(0)
    h=df.meth_diff.abs()
    mod_pval=1-df.ztestpval
    df['val']=h.multiply(mod_pval)
    hh=mod_pval.apply(np.exp)
    exp_val=h.multiply(hh)
    scaled_exp_val=(exp_val.subtract(exp_val.min())).divide(((exp_val.max())-(exp_val.min())))
    
    smooth_exp_val=smoothing(*exp_val) 
    
    scaled_smooth_exp_val=(smooth_exp_val-min(smooth_exp_val))/(max(smooth_exp_val)-min(smooth_exp_val))
    
    df1=pd.concat([df, pd.DataFrame({"exp_val":scaled_exp_val}),pd.DataFrame({"smooth_val":scaled_smooth_exp_val})], axis=1)
    return df1
    
    
def pval_cal_withrep(df):
    filter_col = [col for col in list(df) if col.startswith(('h'))]
    df=df[(df[filter_col]!=0).all(axis=1)]
    
    filter_col1 = [col for col in list(df) if col.startswith(('mc'))]
    filter_col2 = [col for col in list(df) if col.startswith(('h'))]
    mc_read=df[filter_col1]
    total_read=df[filter_col2]
    mc_read.columns=list(total_read.columns.values)
    prop_table=mc_read/total_read
    filter_col1 = [col for col in list(df) if col.startswith(('mc_cont'))]
    filter_col2 = [col for col in list(df) if col.startswith(('mc_case'))]
    prop_names=[]
    DX=[]
    for i in xrange(1,len(filter_col1)+1):
        prop_names.append("meth_cont"+str(i))
        DX.append(1)
    for i in xrange(1,len(filter_col2)+1):
        prop_names.append("meth_case"+str(i))
        DX.append(0)
    prop_table.columns=prop_names
    pval=[]   
    for i in xrange(len(prop_table)):
        props=list(prop_table.ix[i,])
        wgt=list(total_read.ix[i,])
        wgt = [ -x for x in wgt]
            
        tmp_w2 = (1 / (1+np.exp(wgt)))
        weights1= (tmp_w2 - .5) / .5
        data2=pd.DataFrame({'DX':DX,'props':props,'weight':weights1})
        data2.index=range(1,len(data2)+1)
        formula = 'DX ~ props'
        glm_now = ro.r.glm(formula=ro.r(formula), family=ro.r('binomial(link="logit")'), data=data2,weights = weights1)
        res = ro.r.summary(glm_now)
        q=float(np.array(res[7],dtype=float))-float(np.array(res[3],dtype=float))
        d=int(np.array(res[8]))-int(np.array(res[6]))
        if int(np.array(glm_now[18]))==1:
            pval.append(1-(float(np.array(ro.r.pchisq(q, df=d)))))
        else:
            pval.append(1-(1-(float(np.array(ro.r.pchisq(q, df=d))))))
    pval=pd.DataFrame(pval)
    pval.columns=["p_value"]
    df=pd.concat([df,pval],axis=1)
    filter_col3 = [col for col in list(prop_table) if col.startswith(('meth_case'))]
    filter_col4 = [col for col in list(prop_table) if col.startswith(('meth_cont'))]
    
    meth_case_val=prop_table[filter_col3].sum(axis=1)
    meth_cont_val=prop_table[filter_col4].sum(axis=1)
    meth_diff=(meth_case_val/len(filter_col3))-(meth_cont_val/len(filter_col4))
    df['meth_diff']=meth_diff
    df['meth_case']=meth_case_val/len(filter_col3)
    df['meth_cont']=meth_cont_val/len(filter_col4)
    filter_col3 = [col for col in list(df) if col.startswith(('h_case'))]
    filter_col4 = [col for col in list(df) if col.startswith(('h_cont'))]
    h_case_val=df[filter_col3].sum(axis=1)
    h_cont_val=df[filter_col4].sum(axis=1)
    df['h_case']=h_case_val/len(filter_col3)
    df['h_cont']=h_cont_val/len(filter_col4)
    df=df.fillna(0)
    h=df.meth_diff.abs()
    mod_pval=1-df.p_value
    df['val']=h.multiply(mod_pval)
    hh=mod_pval.apply(np.exp)
    exp_val=h.multiply(hh)
    scaled_exp_val=(exp_val-exp_val.min())/(exp_val.max()-exp_val.min())
    
    smooth_exp_val=smoothing(*exp_val) 
    
    scaled_smooth_exp_val=(smooth_exp_val-min(smooth_exp_val))/(max(smooth_exp_val)-min(smooth_exp_val))
    
    df1=pd.concat([df, pd.DataFrame({"exp_val":scaled_exp_val}),pd.DataFrame({"smooth_val":scaled_smooth_exp_val})], axis=1)
    return df1
def smoothing(*a):
    
    import numpy as np
    
    avg_value=[]
    for i in xrange(len(a)):
        if np.sign(i-1)==-1:
            p=0
        else:
            p=a[i-1]
        if i+1==len(a):
            n=0
        else:
            n=a[i+1]
        avg_value.append(np.divide(float(p+a[i]+n),3))
        
    return avg_value
    
def norm_slidingwin_predict_CG(df_file,input_file_path,model_path):
    b=-0.05632
    m=1.89323
    norm_value=[]
    input_file1=input_file_path
    
    df_file1 = pd.read_csv(input_file1,header=None,delimiter=',')

    x=[]
    status=[]
    clf = None
 
    clf = joblib.load(model_path)
    x=np.array(df_file1)

    scaler = preprocessing.StandardScaler().fit(x)
    for i in xrange(len(df_file)-1):
  
            pos_index=i
            
            if (pos_index-5)<0:
                start=pos_index
            else:
                start=pos_index-5
            if (pos_index+5)>=len(df_file):
                stop=len(df_file)-1
            else:
                stop=pos_index+5
            meth_diff=df_file.meth_diff[start:stop]
            value=df_file.smooth_val[start:stop]
            pos_specific=df_file.pos[start:stop]
            
            pos1=df_file.pos[i]
            mod_value=np.ceil(value*10)/10
            sign_win=np.sign(np.median(meth_diff))

            status.append(sign_win)
            val=(abs(pos_specific-pos1)/250.0)
            wght=[]
            for i in val:
                t=min(i,1)
                wght.append(1-t)
            
            
            bins=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            hist,edges = np.histogram(mod_value, bins=bins,weights=wght)
            k=float(sum(hist))
            sum1=hist/k
            norm_value.append(sum1)
    X_test_scaler=scaler.transform(norm_value)
    y_pred=pd.DataFrame(((clf.decision_function(X_test_scaler))), columns=['predicted_values'],dtype='float' )
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    k=pd.concat([df_file.pos[:-1],y_final,status], names=None,axis=1)
    
    return (k)
def norm_slidingwin_predict_nonCG(df_file,input_file_path,model_path):
   
    norm_value=[]
    input_file1=input_file_path
    
    df_file1 = pd.read_csv(input_file1,header=None,delimiter=',')

    x=[]
    status=[]
    clf = None
 
    clf = joblib.load(model_path)
    x=np.array(df_file1)

    scaler = preprocessing.StandardScaler().fit(x)
    for i in xrange(len(df_file)-1):
  
            pos_index=i
            #print pos_index
            if (pos_index-25)<0:
                start=pos_index
            else:
                start=pos_index-25
            if (pos_index+25)>=len(df_file):
                stop=len(df_file)-1
            else:
                stop=pos_index+25
            meth_diff=df_file.meth_diff[start:stop]
            value=df_file.smooth_val[start:stop]
            pos_specific=df_file.pos[start:stop]
            
            pos1=df_file.pos[i]
            mod_value=np.ceil(value*100)/100
            sign_win=np.sign(np.median(meth_diff))

            status.append(sign_win)
            val=(abs(pos_specific-pos1)/10.0)
            wght=[]
            for i in val:
                t=min(i,1)
                wght.append(1-t)
         
            
            bins=np.linspace(0,1,11) 
            hist,edges = np.histogram(mod_value, bins=bins,weights=wght)
            k=float(sum(hist))
            sum1=hist/k
            norm_value.append(sum1)
    X_test_scaler=scaler.transform(norm_value)
    y_pred=pd.DataFrame(((clf.decision_function(X_test_scaler))), columns=['predicted_values'],dtype='float' )
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    k=pd.concat([df_file.pos[:-1],y_pred,status], names=None,axis=1)
    
    return (k)
    
def clustandtrim(k,df1,sc,tr,dis_thres,ncb,score_thres,prn):
    dmr_start=[]
    dmr_stop=[]
    win=False
    no_c=[]
    status=[]
    comb_diff=[]
    mean_diff_case=[]
    mean_diff_cont=[]
    coverage_sample1=[]
    coverage_sample2=[]
    wins=[]
    file_idx=0
    start=0
    stop=0
    for i in xrange(len(k)-1):
        pos1=int(k.pos[i])
        pos2=int(k.pos[i+1])
        label=k.glm_predicted_values[i]
        
        if label>=sc and win==False and (pos2-pos1)<500:
            start=k.pos[i]
            
            win=True
        if win==True and (pos2-pos1)>500:
            
            stop=k.pos[i]
            win=False
        if label<sc and win==True :
            stop=k.pos[i-1]

                       
            win=False  
        if start==stop:
            start=0
            stop=0    
        if start!=0 and stop!=0:    
             dmr_start.append(int(start))
             dmr_stop.append(int(stop))
             p_start=df1[df1.pos==start].index
             p_stop=df1[df1.pos==stop].index
             total_meth_case=df1.meth_case[p_start[0]:p_stop[0]]
             avg_meth_case=total_meth_case.sum()/(len(total_meth_case))
             total_meth_cont=df1.meth_cont[p_start[0]:p_stop[0]]
             avg_meth_cont=total_meth_cont.sum()/(len(total_meth_cont))
             delta_diff=avg_meth_case-avg_meth_cont
             wins.append(np.sign(delta_diff))
             start=0
             stop=0                                           
        
    dmr_start=pd.DataFrame(dmr_start,columns=['dmr_start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['dmr_stop'],dtype='int')     
    final_dmrs=pd.concat([dmr_start,dmr_stop],axis=1)
    dmr_wins=pd.DataFrame(wins,columns=['win_sign'],dtype='int')   
    final_dmrs=pd.concat([dmr_start,dmr_stop,dmr_wins],axis=1)
    dmr_start=[]
    dmr_stop=[]
    m=0
    for i in xrange(len(final_dmrs)-1):
        
       start=final_dmrs.dmr_stop[i]
       stop=final_dmrs.dmr_start[i+1]
       start_indx=df1.pos[df1.pos==start].index[0]
       stop_indx=df1.pos[df1.pos==stop].index[0]
       no_c=((stop_indx-start_indx)-1)
       dmr_dis=stop-start
       k_start=k.pos[k.pos==start].index[0]
       k_stop=k.pos[k.pos==stop].index[0]
       mean_score=np.mean(k.glm_predicted_values[k_start:k_stop+1])
       if dmr_dis>dis_thres:
           new_dmr_start=final_dmrs.dmr_start[i-m]
           new_dmr_stop=final_dmrs.dmr_stop[i]
           m=0
           
           dmr_start.append(int(new_dmr_start))
           dmr_stop.append(int(new_dmr_stop))
       elif final_dmrs.win_sign[i]!=final_dmrs.win_sign[i+1]:
           new_dmr_start=final_dmrs.dmr_start[i-m]
           new_dmr_stop=final_dmrs.dmr_stop[i]
           m=0
           
           dmr_start.append(int(new_dmr_start))
           dmr_stop.append(int(new_dmr_stop))    
       elif dmr_dis<dis_thres and no_c>ncb and mean_score<score_thres:
           new_dmr_start=final_dmrs.dmr_start[i-m]
           new_dmr_stop=final_dmrs.dmr_stop[i]
           m=0
           
           dmr_start.append(int(new_dmr_start))
           dmr_stop.append(int(new_dmr_stop))    
       else:
            m=m+1
    if m!=1:
         dmr_start.append(int(final_dmrs.dmr_start[-1:]))
         dmr_stop.append(int(final_dmrs.dmr_stop[-1:])) 
    dmr_start=pd.DataFrame(dmr_start,columns=['dmr_start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['dmr_stop'],dtype='int')     
    final_dmrs=pd.concat([dmr_start,dmr_stop],axis=1) 
    dmr_start=[]
    dmr_stop=[]
    win=False
    no_c=[]
    status=[]
    comb_diff=[]
    coverage_sample1=[]
    coverage_sample2=[]
    
    file_idx=0
    start=0
    stop=0
    for i in xrange(len(final_dmrs)):  
        start=final_dmrs.dmr_start[i]
        stop=final_dmrs.dmr_stop[i]
        win1=False
        pos=[]
        value=[]
        s=0
        p=0
   
        for j in xrange(file_idx,len(df1)):

            if df1.pos[j]>=start and df1.pos[j]<=stop:
               
               if win1==False:
                  file_idx=j
                  win1=True
               pos.append(df1.pos[j])

               val=df1.exp_val[j]
               
               value.append(val)
            elif df1.pos[j]>stop:
           
                if len(value)>2 :
                   p=len(value)-1
                   while any(t <tr for t in value[s:s+prn]):
                   #while any(t <tr for t in value[s:s+1]):
                     if  s==p:
                         s=p
                         break
                     else:
                         s=s+1
                
                   cg_start=pos[s]
                
                   while any(t <tr for t in value[p-(prn-1):p+1]):
                   #while any(t <tr for t in value[p:p+1]):    
                       if p==(prn-1):
                       #if p==1:
                           p=len(value)-1
                           break
                       else:
                          p=p-1
                
                   cg_stop=pos[p]
                   
                   if (cg_stop-cg_start)>=10:
                       p_start=df1[df1.pos==cg_start].index
                       p_stop=df1[df1.pos==cg_stop].index
                       dif=df1.meth_diff[p_start[0]:p_stop[0]]
                       
                       coverage_cont=np.median(df1.h_cont[p_start[0]:p_stop[0]])
                       coverage_case=np.median(df1.h_case[p_start[0]:p_stop[0]])
                       total_meth_case=df1.meth_case[p_start[0]:p_stop[0]]
                       avg_meth_case=total_meth_case.sum()/(len(total_meth_case))
                       total_meth_cont=df1.meth_cont[p_start[0]:p_stop[0]]
                       avg_meth_cont=total_meth_cont.sum()/(len(total_meth_cont))
                       delta_diff=avg_meth_case-avg_meth_cont
                       dif=dif.abs()
                      
                       start_indx=k.pos[k.pos==cg_start].index[0]
                       stop_indx=k.pos[k.pos==cg_stop].index[0]
                       win_sign=k.win_sign[start_indx:stop_indx+1]
                       no_c.append(len(win_sign))

                       if np.sign(delta_diff)==-1:
                           status_win="hypo"
                       if np.sign(delta_diff)==1:
                           status_win="hyper"
                       if np.sign(delta_diff)==0:
                           status_win="hemi"    
                       status.append(status_win)
                       dmr_start.append(cg_start)
                       dmr_stop.append(cg_stop)
                       comb_diff.append(delta_diff) 
                       mean_diff_case.append(avg_meth_case)
                       mean_diff_cont.append(avg_meth_cont)
                       coverage_sample1.append(coverage_cont)
                       coverage_sample2.append(coverage_case)
                break 
    mean_diff_case=pd.DataFrame(mean_diff_case,columns=['mean_Meth1'],dtype='float')
    mean_diff_cont=pd.DataFrame(mean_diff_cont,columns=['mean_Meth2'],dtype='float')              
    comb_diff=pd.DataFrame(comb_diff,columns=['delta'],dtype='float')
    dmr_start=pd.DataFrame(dmr_start,columns=['start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['end'],dtype='int')
    status_dmr=pd.DataFrame(status,columns=['status'],dtype='str')
    number_of_Cs=pd.DataFrame(no_c,columns=['numCG'],dtype='int')
    avg_number_of_Cs_sample1=pd.DataFrame(coverage_sample1,columns=['avg_coverage_of_Cs_sample1'],dtype='int')
    avg_number_of_Cs_sample2=pd.DataFrame(coverage_sample2,columns=['avg_coverage_of_Cs_sample2'],dtype='int')
    length=pd.DataFrame(dmr_stop.end-dmr_start.start,columns=['length'],dtype='int')
    final_dmrs=pd.concat([dmr_start,dmr_stop,status_dmr,number_of_Cs,mean_diff_case,mean_diff_cont,comb_diff,avg_number_of_Cs_sample1,avg_number_of_Cs_sample2,length],axis=1)
    
    final_dmrs=final_dmrs.reset_index(drop=True)
    return (final_dmrs)
    
def filterdmr(dmrs,minlen,mc,d):
    final_dmrs=dmrs.loc[(dmrs.numCG >= mc) & (abs(dmrs.delta)>= d) & (abs(dmrs.length)>= minlen)]
    final_dmrs=final_dmrs.reset_index(drop=True)
    return (final_dmrs)