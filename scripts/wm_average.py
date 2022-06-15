# Script for importing the MDTB data set from super_cerebellum to general format.
import pandas as pd
import numpy as np
import nibabel as nb
import SUITPy as suit
import matplotlib.pyplot as plt
import seaborn as sb
import atlas_map as am 
import surfAnalysisPy as surf

base_dir = '/Volumes/diedrichsen_data$/data/Cerebellum/CerebellumWorkingMemory'

roi_dir = base_dir + '/rois'
dataC_dir = base_dir + '/WM_maps/cortex'

def avrg_profiles():
    # Make the atlas object
    surf_roi =[] 
    hem_name = ['cortex_left','cortex_right']
    avrg_data = []
    # Get the parcelation 
    for i,h in enumerate(['L','R']): 
        gifti = roi_dir + f'/rois_hand.{h}.label.gii'
        LG = nb.load(gifti)
        LG.labeltable.print_summary()
        surf_roi.append(am.AtlasSurfaceParcel(hem_name[i],gifti))

        datag1 = dataC_dir + f'/wgroup.encode-rest.{h}.func.gii'
        datag2 = dataC_dir + f'/wgroup.retriev-rest.{h}.func.gii'
        G1 = nb.load(datag1)
        G2 = nb.load(datag2)
        Data = np.c_[np.c_[G1.agg_data()],np.c_[G2.agg_data()]]
        avrg_data.append(surf_roi[i].agg_data(Data.T))
    pass
    data = np.c_[avrg_data[0],avrg_data[1]]
    # Make data structure 
    T=pd.DataFrame(data)
    rois = np.array(['LVis','LTrans','LIP','LM1','LvPMD','LrPMD','RVis','RTrans','RIP','RM1','RvPMD','RrPMD'])
    T['load'] = [2,4,6,2,4,6,2,4,6,2,4,6]
    T['dir'] = ['B','B','B','F','F','F','B','B','B','F','F','F']
    T['phase'] = ['E','E','E','E','E','E','R','R','R','R','R','R']
    D= pd.melt(T,id_vars=['load','dir','phase'])
    D['roi'] = rois[D.variable.to_numpy(dtype=int)]
    for i in range(12):
        plt.subplot(4,3,i+1)
        if i==0:
            sb.lineplot(data=D[D.variable==i],x='load',y='value',hue='dir',style='phase')
        else:
            sb.lineplot(data=D[D.variable==i],x='load',y='value',hue='dir',style='phase',legend=False)
        plt.ylim([0.1,0.60])
        plt.xlabel(None)
        plt.title(rois[i])
    pass



if __name__ == "__main__":
    avrg_profiles()
    pass


    # T= pd.read_csv(data_dir + '/participants.tsv',delimiter='\t')
    # for s in T.participant_id:
    #     pass