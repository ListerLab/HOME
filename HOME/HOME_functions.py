
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 19:06:05 2017

@author: akanksha
"""
import subprocess
import pandas as pd 
from itertools import groupby
import numpy as np

from sklearn import preprocessing
#from sklearn.externals import joblib
import statsmodels.stats.proportion as sm
def fun_win(val):
     t=1-(min(val,1))
     return t 

def fill_na(df_file):
    filter_col = [col for col in list(df_file) if col.startswith(('mc'))]
    filter_col1 = [col for col in list(df_file) if col.startswith(('h'))]
    df=df_file[filter_col]
    df=df.bfill(axis=1).ffill(axis=1)
    df_file[filter_col]
    df_file[filter_col]=df
    df_file[filter_col]=df_file[filter_col].astype(int)
    
    df=df_file[filter_col1]
    df=df.bfill(axis=1).ffill(axis=1)
    df_file[filter_col1]=df
    df_file[filter_col1]=df_file[filter_col1].astype(int)
    return df_file

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
        df=df[(df[filter_col]> 0).all(axis=1)]
        df=df.reset_index(drop=True)
   elif classes=="CHN" or classes=="CHG" or classes=="CHH" or classes=="CNN":

        filter_col = [col for col in list(df) if col.startswith(('chr','pos','mc','h'))]
        df=df[filter_col]
        df=df.reset_index(drop=True)
        filter_col = [col for col in list(df) if col.startswith(('h'))]
        df=df[(df[filter_col]> 0).all(axis=1)]
        df=df.reset_index(drop=True)        
      
   return df
def process_frame_withR(file1):
    com="Rscript ./scripts/HOME_R.R" + " "+file1
    subprocess.call(com, shell=True)
    
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
    df1=df1.fillna(0)  
    return df1
def chunker1(seq, sizes):
 idx=0
 for s in sizes:

    k=seq[idx:idx+s]
    idx=idx+s
 
    yield k
    
def pval_format_withrep(df_path):
    df=pd.read_table(df_path,header=0)
    
    filter_col = [col for col in list(df) if col.startswith(('h'))]
    df=df[(df[filter_col]!=0).all(axis=1)]
#    
    filter_col1 = [col for col in list(df) if col.startswith(('mc'))]
    filter_col2 = [col for col in list(df) if col.startswith(('h'))]
    mc_read=df[filter_col1]
    total_read=df[filter_col2]
    mc_read.columns=list(total_read.columns.values)
    prop_table=mc_read/total_read
    del(mc_read)
    filter_col1 = [col for col in list(df) if col.startswith(('mc_cont'))]
    filter_col2 = [col for col in list(df) if col.startswith(('mc_case'))]
    prop_names=[]
   
    for i in xrange(1,len(filter_col1)+1):
        prop_names.append("meth_cont"+str(i))
        
    for i in xrange(1,len(filter_col2)+1):
        prop_names.append("meth_case"+str(i))
        
    prop_table.columns=prop_names
    filter_col3 = [col for col in list(prop_table) if col.startswith(('meth_case'))]
    filter_col4 = [col for col in list(prop_table) if col.startswith(('meth_cont'))]
    
    meth_case_val=prop_table[filter_col3].sum(axis=1)
    meth_cont_val=prop_table[filter_col4].sum(axis=1)
    df['meth_case']=meth_case_val/len(filter_col3)
    df['meth_cont']=meth_cont_val/len(filter_col4)
    meth_diff=df.meth_case-df.meth_cont
    df['meth_diff']=meth_diff

    filter_col1 = [col for col in list(df) if col.startswith(('h_case'))]
    filter_col2= [col for col in list(df) if col.startswith(('h_cont'))]
    h_case_val=df[filter_col1].sum(axis=1)
    h_cont_val=df[filter_col2].sum(axis=1)
    df['h_case']=h_case_val/len(filter_col1)
    df['h_cont']=h_cont_val/len(filter_col2)
    df=df.fillna(0)
    h=df.meth_diff.abs()
    mod_pval=1-df.p_value
    df['val']=h.multiply(mod_pval)
    hh=mod_pval.apply(np.exp)
    df['exp_val1']=h.multiply(hh)
    df['exp_val']=(df.exp_val1-df.exp_val1.min())/(df.exp_val1.max()-df.exp_val1.min())
    df=df.fillna(0) 
    return df
    
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
def chunker(seq, size):
  
   for pos in xrange(0, len(seq), size):
    
    start_df=max(0,pos-25)
    start=max(0,pos)
    stop_df=min(len(seq),pos+size+25)
    stop=min(len(seq),((pos+size)-1))
    k=seq[start_df:stop_df]
    
  
    yield (k,start,stop)    
def norm_slidingwin_predict_CG(df_file,input_file_path,model_path):
    b=-0.123176135253
    m=1.95258046977
#    b=-0.05632
#    m=1.89323
    norm_value=[]
    #input_file1=input_file_path
    x=input_file_path
    #df_file1 = pd.read_csv(input_file1,header=None,delimiter=',')

    
    status=[]
    #clf = None
    W=np.load(model_path+"W.npy")
    bb=np.load(model_path+"b.npy")
    #clf = joblib.load(model_path)
   # x=np.array(df_file1)

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
            
            wght=map(fun_win, val)
            
            bins=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            hist,edges = np.histogram(mod_value, bins=bins,weights=wght)
            k=float(sum(hist))
            sum1=hist/k
            norm_value.append(sum1)
    X_test_scaler=scaler.transform(norm_value)
    y_pred_int=np.dot(X_test_scaler,np.transpose(W))+bb
    y_pred=pd.DataFrame(y_pred_int, columns=['predicted_values'],dtype='float' )
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    k=pd.concat([df_file.pos[:-1],y_final,status], names=None,axis=1)
    
    return (k)
def norm_slidingwin_predict_nonCG(df_file,input_file_path,model_path):
     
#    b=0.57559 
#    m=2.01748
    b=0.21610910497
    m=1.01167212729
    norm_value=[]
    #input_file1=input_file_path
    x=input_file_path
   # df_file1 = pd.read_csv(input_file1,header=None,delimiter=',')
    
    
    status=[]
    #clf = None
    W=np.load(model_path+"W.npy")
    bb=np.load(model_path+"b.npy") 
    #clf = joblib.load(model_path)
    #x=np.array(df_file1)
    
    scaler = preprocessing.StandardScaler().fit(x)
    
    for i in df_file[0].index:
        if i>=df_file[1] and i<=df_file[2]: 
            pos_index=i
            start=max(0,pos_index-25)
            stop=min(df_file[2],pos_index+25)
            
            meth_diff=df_file[0].loc[start:stop].meth_diff
            value=df_file[0].loc[start:stop].smooth_val
            pos_specific=df_file[0].loc[start:stop].pos
            
            pos1=df_file[0].pos[i]
            mod_value=np.ceil(value*100)/100
            sign_win=np.sign(np.median(meth_diff))
    
            status.append(sign_win)
            val=(abs(pos_specific-pos1)/10.0)
            wght=map(fun_win, val)
         
            bins=np.linspace(0,1,11) 
            hist,edges = np.histogram(mod_value, bins=bins,weights=wght)
            k=float(sum(hist))
            sum1=hist/k
            norm_value.append(sum1)
        if i>df_file[2]:  
             break
    #print norm_value
    X_test_scaler=scaler.transform(norm_value)
    y_pred=pd.DataFrame(np.dot(X_test_scaler,np.transpose(W))+bb) 
    y_pred.columns=['predicted_values']
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    df=df_file[0].pos.ix[df_file[1]:df_file[2]].reset_index(drop=True)
    k=pd.concat([df,y_final,status], names=None,axis=1) 

    return (k)

    
def clustandtrim_CG(k,df1,sc,tr,dis_thres,ncb,prn,len_cutoff):
    score_thres=0.02
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
    if len(final_dmrs)>1:
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
                       
                       if (cg_stop-cg_start)>=len_cutoff:
                           p_start=df1[df1.pos==cg_start].index
                           p_stop=df1[df1.pos==cg_stop].index
                           dif=df1.meth_diff[p_start[0]:p_stop[0]]
                           
                           coverage_cont=np.median(df1.h_cont[p_start[0]:p_stop[0]])
                           coverage_case=np.median(df1.h_case[p_start[0]:p_stop[0]])
                           total_meth_case=df1.meth_case[p_start[0]:p_stop[0]]
                           avg_meth_case=total_meth_case.sum()/(len(total_meth_case))
                           total_meth_cont=df1.meth_cont[p_start[0]:p_stop[0]]
                           avg_meth_cont=total_meth_cont.sum()/(len(total_meth_cont))
                           delta_diff=avg_meth_cont-avg_meth_case
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
        mean_diff_case=pd.DataFrame(mean_diff_case,columns=['mean_Meth2'],dtype='float')
        mean_diff_cont=pd.DataFrame(mean_diff_cont,columns=['mean_Meth1'],dtype='float')              
        comb_diff=pd.DataFrame(comb_diff,columns=['delta'],dtype='float')
        dmr_start=pd.DataFrame(dmr_start,columns=['start'],dtype='int')
        dmr_stop=pd.DataFrame(dmr_stop,columns=['end'],dtype='int')
        status_dmr=pd.DataFrame(status,columns=['status'],dtype='str')
        number_of_Cs=pd.DataFrame(no_c,columns=['numC'],dtype='int')
        avg_number_of_Cs_sample1=pd.DataFrame(coverage_sample1,columns=['avg_coverage1'],dtype='int')
        avg_number_of_Cs_sample2=pd.DataFrame(coverage_sample2,columns=['avg_coverage2'],dtype='int')
        length=pd.DataFrame(dmr_stop.end-dmr_start.start,columns=['len'],dtype='int')
        final_dmrs=pd.concat([dmr_start,dmr_stop,status_dmr,number_of_Cs,mean_diff_cont,mean_diff_case,comb_diff,avg_number_of_Cs_sample1,avg_number_of_Cs_sample2,length],axis=1)
        final_dmrs=final_dmrs.loc[(final_dmrs.status != 'hemi')]
        final_dmrs=final_dmrs.reset_index(drop=True)
        return (final_dmrs)
    else:
        final_dmrs=[]
        return final_dmrs
def clustandtrim_nonCG1(k,sc):
    k.reset_index(drop=True,inplace=True)
    dmr_start=[]
    dmr_stop=[]
    win=False
    wins=[]
    num_c=[]
 
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
  
        if start!=0 and stop!=0:    
             dmr_start.append(int(start))   
             dmr_stop.append(int(stop))
             
             p_start=k[k.pos==start].index
             p_stop=k[k.pos==stop].index
             total_meth_case=k.meth_case[p_start[0]:p_stop[0]+1]
             avg_meth_case=total_meth_case.sum()/(len(total_meth_case))
             total_meth_cont=k.meth_cont[p_start[0]:p_stop[0]+1]
             avg_meth_cont=total_meth_cont.sum()/(len(total_meth_cont))
             delta_diff=avg_meth_case-avg_meth_cont
             num_c.append(len(total_meth_case))
             wins.append(np.sign(delta_diff))
             
             start=0
             stop=0                                           
        
    dmr_start=pd.DataFrame(dmr_start,columns=['dmr_start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['dmr_stop'],dtype='int')     
    final_dmrs=pd.concat([dmr_start,dmr_stop],axis=1)
    dmr_wins=pd.DataFrame(wins,columns=['win_sign'],dtype='int')
    numb_c=pd.DataFrame(num_c,columns=['numc'],dtype='int') 
    final_dmrs=pd.concat([dmr_start,dmr_stop,dmr_wins,numb_c],axis=1)
    if len (final_dmrs)>1:
        return (final_dmrs)
   
    else:
        final_dmrs=pd.DataFrame()
        return final_dmrs
    
    
def splitlist(k,df,npp,dis_thres):
    h=df.dmr_start[1:]
    h.reset_index(drop=True,inplace=True)
    hh=df.dmr_stop[0:-1]
    l=h-hh
    l.loc[-1]=0
    l.index = l.index + 1
    l.sort_index(inplace=True)
    df["dis"]=l
    
    chunks= int(len(df)/npp)
    chunks_list=[chunks]*npp 
    chunks_list=np.cumsum(chunks_list)
    chunks_list=chunks_list-chunks_list[0]
   
    mod_chunks_list=[]
    for i in chunks_list:
        
        c=0
        if df.dis[i]>dis_thres or df.dis[i]==0 :
            mod_chunks_list.append(i)
        elif df.dis[i]<dis_thres:
            
            while (df.dis[(i+c)]<dis_thres) and (i+c<=len(df)-2):
                
                c=c+1
                
            mod_chunks_list.append(i+c)
                 
    if    mod_chunks_list[-1]!=len(df)-1:
        mod_chunks_list.append(len(df)-1)
    
    mod_chunks_list=[ key for key,_ in groupby(mod_chunks_list)]
    
    for i in xrange(len(mod_chunks_list)-1):
       
       if i== len(mod_chunks_list)-2:
        df_split=df.ix[mod_chunks_list[i]:mod_chunks_list[i+1]]
       else:
           df_split=df.ix[mod_chunks_list[i]:mod_chunks_list[i+1]-1]
       df_split.reset_index(drop=True,inplace=True)
       
       k_split=k.ix[k.pos[k.pos==df_split.dmr_start[0]].index[0]:k.pos[k.pos==df_split.dmr_stop[len(df_split)-1]].index[0]]
       k_split.reset_index(drop=True,inplace=True)
       yield  df_split,k_split
       
def clustandtrim_nonCG2(k,final_dmrs,dis_thres,ncb,len_cutoff):
    final_dmrs.reset_index(drop=True,inplace=True)
 
    dmr_start=[]
    dmr_stop=[]
    
    m=1
    
    for i in xrange(1,len(final_dmrs)):
        
       start=final_dmrs.dmr_stop[i-1]
       stop=final_dmrs.dmr_start[i]
       numc_prev=final_dmrs.numc[i-1]
       numc_current=final_dmrs.numc[i]
       len_prev=final_dmrs.dmr_stop[i-1]-final_dmrs.dmr_start[i-1]
       len_current=final_dmrs.dmr_stop[i]-final_dmrs.dmr_start[i]
       dmr_dis=stop-start
    
       
       if dmr_dis>dis_thres:
           new_dmr_start=final_dmrs.dmr_start[i-m]
           new_dmr_stop=final_dmrs.dmr_stop[i-1]
           m=1
           dmr_start.append(int(new_dmr_start))
           dmr_stop.append(int(new_dmr_stop))
       elif final_dmrs.win_sign[i]!=final_dmrs.win_sign[i-1] and numc_current>5 and numc_prev>5:
         if  len_prev>=(len_current/2) and len_current>=(len_prev/2):
           new_dmr_start=final_dmrs.dmr_start[i-m]
           new_dmr_stop=final_dmrs.dmr_stop[i-1]
           m=1
           
           dmr_start.append(int(new_dmr_start))
           dmr_stop.append(int(new_dmr_stop))    

           
       else:
           
            m=m+1
    if m!=1:
           new_dmr_start=final_dmrs.dmr_start[i-(m-1)]
           new_dmr_stop=final_dmrs.dmr_stop[i]
           dmr_start.append(int(new_dmr_start))
           dmr_stop.append(int(new_dmr_stop))
    if i==(len(final_dmrs)-1) and m==1:
          
           dmr_start.append(int(final_dmrs.dmr_start[-1:]))
           dmr_stop.append(int(final_dmrs.dmr_stop[-1:]))        
       
         
    dmr_start=pd.DataFrame(dmr_start,columns=['dmr_start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['dmr_stop'],dtype='int')     
    final_dmrs=pd.concat([dmr_start,dmr_stop],axis=1)

    dmr_start=[]
    dmr_stop=[]
    
    no_c=[]
    status=[]
    comb_diff=[]
    coverage_sample1=[]
    coverage_sample2=[]
    mean_diff_case=[]
    mean_diff_cont=[]
    
    for i in xrange(len(final_dmrs)):
        cg_start=final_dmrs.dmr_start[i]
        cg_stop=final_dmrs.dmr_stop[i]
        
        if (cg_stop-cg_start)>=len_cutoff:
                   
                   p_start=k[k.pos==cg_start].index
                   p_stop=k[k.pos==cg_stop].index
                   dif=k.meth_diff[p_start[0]:p_stop[0]+1]
                   
                   coverage_cont=np.median(k.h_cont[p_start[0]:p_stop[0]+1])
                   coverage_case=np.median(k.h_case[p_start[0]:p_stop[0]+1])
                   total_meth_case=k.meth_case[p_start[0]:p_stop[0]+1]
                   avg_meth_case=total_meth_case.sum()/(len(total_meth_case))
                   total_meth_cont=k.meth_cont[p_start[0]:p_stop[0]+1]
                   avg_meth_cont=total_meth_cont.sum()/(len(total_meth_cont))
                   delta_diff=avg_meth_cont-avg_meth_case
                   dif=dif.abs()
                  
                   win_sign=k.win_sign[p_start[0]:p_stop[0]+1]
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
                   
    mean_diff_case=pd.DataFrame(mean_diff_case,columns=['mean_Meth2'],dtype='float')
    mean_diff_cont=pd.DataFrame(mean_diff_cont,columns=['mean_Meth1'],dtype='float')              
    comb_diff=pd.DataFrame(comb_diff,columns=['delta'],dtype='float')
    dmr_start=pd.DataFrame(dmr_start,columns=['start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['end'],dtype='int')
    status_dmr=pd.DataFrame(status,columns=['status'],dtype='str')
    number_of_Cs=pd.DataFrame(no_c,columns=['numC'],dtype='int')
    avg_number_of_Cs_sample1=pd.DataFrame(coverage_sample1,columns=['avg_coverage1'],dtype='int')
    avg_number_of_Cs_sample2=pd.DataFrame(coverage_sample2,columns=['avg_coverage2'],dtype='int')
    length=pd.DataFrame(dmr_stop.end-dmr_start.start,columns=['len'],dtype='int')
    final_dmrs=pd.concat([dmr_start,dmr_stop,status_dmr,number_of_Cs,mean_diff_cont,mean_diff_case,comb_diff,avg_number_of_Cs_sample1,avg_number_of_Cs_sample2,length],axis=1)
    final_dmrs=final_dmrs.loc[(final_dmrs.status != 'hemi')]
    final_dmrs=final_dmrs.reset_index(drop=True)
    if len (final_dmrs)>1:
        return (final_dmrs)
   
    else:
        final_dmrs=pd.DataFrame()
        return final_dmrs
   
def norm_slidingwin_predict_nonCG_withoutchunk(df_file,input_file_path,model_path):
     
    b=0.21610910497
    m=1.01167212729
    norm_value=[]
    x=input_file_path
   
    status=[]
    W=np.load(model_path+"W.npy")
    bb=np.load(model_path+"b.npy") 
#    clf = None
# 
#    clf = joblib.load(model_path)
    

    scaler = preprocessing.StandardScaler().fit(x)
    for i in xrange(len(df_file)-1):
  
            pos_index=i
            
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
            wght=map(fun_win, val)
            
            bins=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
            hist,edges = np.histogram(mod_value, bins=bins,weights=wght)
            k=float(sum(hist))
            sum1=hist/k
            norm_value.append(sum1)
    X_test_scaler=scaler.transform(norm_value)
    y_pred_int=np.dot(X_test_scaler,np.transpose(W))+bb 
    y_pred=pd.DataFrame(y_pred_int, columns=['predicted_values'],dtype='float' )
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    k=pd.concat([df_file.pos[:-1],y_final,status], names=None,axis=1)
    
    return (k)
    
def filterdmr(dmrs,minlen,mc,d):
    final_dmrs=dmrs.loc[(dmrs.numC >= mc) & (abs(dmrs.delta)>= d) & (abs(dmrs.len)>= minlen)]
    final_dmrs=final_dmrs.reset_index(drop=True)
    return (final_dmrs)

def filterdmr_nonCG(dmrs,minlen,mc,d):
    final_dmrs=dmrs.loc[(dmrs.numC >= mc) & (abs(dmrs.len)>= minlen)]
    final_dmrs=final_dmrs.reset_index(drop=True)
    return (final_dmrs)