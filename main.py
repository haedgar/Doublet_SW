
# <font color='Red'> $\;$ Reservoir model: a geothermal doublet in a single well</font>
# <font color='Red'> $\;$ Based on original script provided in First open-DARTS workshop 2023</font>
# <font color='blue'> $\;$ https://darts.citg.tudelft.nl/workshops/</font>
# <font color='blue'> $\;$ https://gitlab.com/open-darts/darts-workshop/-/blob/main/2.Geothermal_model.ipynb</font>

# DARTS version 1.0.2

# Edgar Hernandez

import os 

from Model import Model

#vertical permeability sensitivity at different flow rates
# remember to change the parent directory to your local

kvkh = [0.2,0.5,1]      # vertical to horizontal permeability ratio [dimensionless]
flowrates = [30,45,60]  #[m3/hr] injection flow rate
k = 1
for l in kvkh:
    for j in flowrates:
        
        directory = "Caso_"+str(k)
        parent_dir = "C:/Users/HERNANDE/OneDrive - VITO/Documents/Single Well/Concept 2/Real case darts sim/Final sim"
        m = Model(l)
        
        m.init()
        
        # output initial conditions
        out_vtk= 'case_'+str(l) +'_'+ str(j)
        m.export_pro_vtk(out_vtk)
        m.run_python(3650)  # TIME FOR EQUILIBRATION
        #m.run_python(1000)  # TIME FOR EQUILIBRATION
        
        #Temperatures = [38.7+273.15, 37.8+273.15,37.6+273.15,37.2+273.15] # This can be used for changing temeprature overtime
        Dtimes = [365*5,3650,2*3650-365*5-3650]               # delta production times [days]
        #Dtimes = [10]               # delta production times [days]
        
        i = 0
        for t in range(len(Dtimes)):
        #    # run and output every 10 years (30 in total)
            #m.run_python(10*365, restart_dt=0.5) # restart_dt >> time step size after equilibration
            m.set_boundary_conditions(T=40+273.15,Q_INJ = j,BHP = 260)
            m.run_python(Dtimes[i], restart_dt=0.5) # restart_dt >> time step size after equilibratio
            
            m.export_pro_vtk(out_vtk)
            i = i + 1
        
        # print timers and statistics for the run
        m.print_timers()
        m.print_stat()
        
        path = os.path.join(parent_dir, directory)
        
        if os.path.isdir(path) == False:
            os.mkdir(path)
        
        file_name = 'well_data_base_kv_kh_'+str(l)+'_'+str(j)+'.xlsx'
        address = parent_dir+'/'+directory+'/'+file_name
        
        import pandas as pd
        # output well information to Excel file
        td = pd.DataFrame.from_dict(m.physics.engine.time_data)
        
        #writer = pd.ExcelWriter('well_data_base_kv_kh_1_perm_av.xlsx')
        writer = pd.ExcelWriter(address)
        td.to_excel(writer, 'Sheet1')
        writer.save()
    k+=1