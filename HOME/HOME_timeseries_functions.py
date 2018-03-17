import pandas as pd 

import numpy as np
import subprocess
from sklearn import preprocessing
#from sklearn.externals import joblib
import statsmodels.stats.proportion as sm

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
        
        filter_col = [col for col in list(df) if col.startswith(('chr','pos','mc','h'))]
        df=df[filter_col]
        df=df.reset_index(drop=True)        
      
   return df
def pval_cal_withoutrep(df):
    df=df.rename(columns = {'mc2_rep1':'mc_case','mc1_rep1':'mc_cont','h2_rep1':'h_case','h1_rep1':'h_cont'})
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
    
def process_frame_withR(file1):
    com="Rscript ./scripts/HOME_R_time.R" + " "+file1
    subprocess.call(com, shell=True)
    
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
    filter_col1 = [col for col in list(df) if col.startswith(('h1'))]
    filter_col2 = [col for col in list(df) if col.startswith(('h2'))]
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
    df['exp_val']=h.multiply(hh)
    df['scaled_exp_val']=(df.exp_val-df.exp_val.min())/(df.exp_val.max()-df.exp_val.min())
    
    return df

def chunker1(seq, sizes):
 idx=0
 for s in sizes:

    k=seq[idx:idx+s]
    idx=idx+s
 
    yield k
def chunker(seq, size):
  
   for pos in xrange(0, len(seq), size):
    
    start_df=max(0,pos-25)
    start=max(0,pos)
    stop_df=min(len(seq),pos+size+25)
    stop=min(len(seq),((pos+size)-1))
    k=seq[start_df:stop_df]
    
  
    yield (k,start,stop)        
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
    delta=[]
    x=[]
    status=[]
    W=np.load(model_path+"W.npy")
    bb=np.load(model_path+"b.npy")
#    clf = None
# 
#    clf = joblib.load(model_path)
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
            mc_diff=df_file.meth_diff[i]
            delta.append(mc_diff)
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
    y_pred=pd.DataFrame(np.dot(X_test_scaler,np.transpose(W))+bb)
    
    y_pred.columns=['predicted_values']
    
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    delta=pd.DataFrame(delta,columns=['delta'],dtype='float')
    k=pd.concat([df_file.pos[:-1],y_final,delta,status], names=None,axis=1)       
   
    
    return (k) 
def norm_slidingwin_predict_nonCG_withoutchunk(df_file,input_file_path,model_path):
     
    b=0.57559 
    m=2.01748
    norm_value=[]
    input_file1=input_file_path
    
    df_file1 = pd.read_csv(input_file1,header=None,delimiter=',')

    x=[]
    status=[]
    #clf = None
    delta=[]
    W=np.load(model_path+"W.npy")
    bb=np.load(model_path+"b.npy")
    #clf = joblib.load(model_path)
    x=np.array(df_file1)

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
            delta.append(df_file.meth_diff[i])
            pos1=df_file.pos[i]
            mod_value=np.ceil(value*100)/100
            sign_win=np.sign(np.median(meth_diff))

            status.append(sign_win)
            val=(abs(pos_specific-pos1)/10.0)
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
    y_pred_int=np.dot(X_test_scaler,np.transpose(W))+bb
    y_pred=pd.DataFrame(y_pred_int, columns=['predicted_values'],dtype='float' )
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    
    delta=pd.DataFrame(delta,columns=['delta'],dtype='float')
    
    k=pd.concat([df_file.pos[:-1],y_final,delta,status], names=None,axis=1)
    
    return (k)
 
def norm_slidingwin_predict_nonCG(df_file,input_file_path,model_path):
     
    b=0.57559 
    m=2.01748
    norm_value=[]
    input_file1=input_file_path
    
    df_file1 = pd.read_csv(input_file1,header=None,delimiter=',')
    delta=[]
    x=[]
    status=[]
    #clf = None
    W=np.load(model_path+"W.npy")
    bb=np.load(model_path+"b.npy") 
    #clf = joblib.load(model_path)
    x=np.array(df_file1)
    
    scaler = preprocessing.StandardScaler().fit(x)
    
    for i in df_file[0].index:
        if i>=df_file[1] and i<=df_file[2]: 
            pos_index=i
            start=max(0,pos_index-25)
            stop=min(df_file[2],pos_index+25)
            
            meth_diff=df_file[0].loc[start:stop].meth_diff
            value=df_file[0].loc[start:stop].smooth_val
            pos_specific=df_file[0].loc[start:stop].pos
            mc_diff=df_file[0].meth_diff[i]
            delta.append(mc_diff)
            pos1=df_file[0].pos[i]
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
        if i>df_file[2]:  
             break
    #print norm_value
    X_test_scaler=scaler.transform(norm_value)
    y_pred=pd.DataFrame(np.dot(X_test_scaler,np.transpose(W))+bb)
    y_pred.columns=['predicted_values']
    y_final= np.exp((b + m*y_pred)) / (1 + np.exp((b + m*y_pred))) 
    y_final.columns=['glm_predicted_values']
    status=pd.DataFrame(status,columns=['win_sign'],dtype='float')
    delta=pd.DataFrame(delta,columns=['delta'],dtype='float')
    df=df_file[0].pos.ix[df_file[1]:df_file[2]].reset_index(drop=True)
    k=pd.concat([df,y_final,delta,status], names=None,axis=1)  

    return (k)
def clustandtrim(k,sc,minlen,mc):
    dmr_start=[]
    dmr_stop=[]
    win=False
    no_c=[]
    start=0
    stop=0
    for i in xrange(len(k)-1):
        pos1=int(k.pos[i])
        pos2=int(k.pos[i+1])
        label=k.glm_predicted_values[i]
        #print label
        if label>=sc and win==False and (pos2-pos1)<500:
            start=k.pos[i]
            
            
            win=True            
        if win==True and (pos2-pos1)>500:
            
            stop=k.pos[i]

            
            
        if label<sc and win==True :
            stop=k.pos[i-1]

            win=False
        

        if start!=0 and stop!=0:    
            dmr_start.append(start)
            dmr_stop.append(stop)
            start_indx=k.pos[k.pos==start].index[0]
            stop_indx=k.pos[k.pos==stop].index[0]
            win_len=k.pos[start_indx:stop_indx+1]
            no_c.append(len(win_len))
            win=False
            start=0
            stop=0
    dmr_start=pd.DataFrame(dmr_start,columns=['start'],dtype='int')
    dmr_stop=pd.DataFrame(dmr_stop,columns=['end'],dtype='int')
    length=pd.DataFrame(dmr_stop.end-dmr_start.start,columns=['len'],dtype='int')
    number_of_Cs=pd.DataFrame(no_c,columns=['numC'],dtype='int')
    
    final_dmrs=pd.concat([dmr_start,dmr_stop,number_of_Cs,length],axis=1)
    final_dmrs=final_dmrs[final_dmrs.len>=minlen]
    final_dmrs=final_dmrs[final_dmrs.numC>=mc]
    final_dmrs=final_dmrs.reset_index(drop=True)
    return (final_dmrs)    