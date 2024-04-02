#DARTS version 1.0.2

# <font color='Red'> $\;$ Reservoir model: a geothermal doublet in a single well</font>
# <font color='Red'> $\;$ Based on original script provided in First open-DARTS workshop 2023</font>
# <font color='blue'> $\;$ https://darts.citg.tudelft.nl/workshops/</font>
# <font color='blue'> $\;$ https://gitlab.com/open-darts/darts-workshop/-/blob/main/2.Geothermal_model.ipynb</font>
# Edgar Hernandez


from darts.models.physics.geothermal import Geothermal                 # define number of component
from darts.models.darts_model import DartsModel
from darts.models.physics.iapws.iapws_property_vec import _Backward1_T_Ph_vec  # transformation from T to entalpy
from darts.tools.keyword_file_tools import load_single_keyword          # TO RECOGNICE SOME KEYWORDS
from darts.engines import redirect_darts_output
import numpy as np
from darts.tools.keyword_file_tools import load_single_keyword          # TO RECOGNICE SOME KEYWORDS
redirect_darts_output('run_geothermal.log')
from iapws.iapws97 import _Region1

## Check property container

# Reading por and perm and editing order of input in the case they are coming from Sgems geostatistic package where data is ordered from the bottom to the top of the model
filename3 = 'por.txt'
filename4 = 'perm_h.txt'
poro_INV = load_single_keyword(filename3, 'Por')  # order of permeability from Sgem is inverse
perm_INV = load_single_keyword(filename4, 'KH_P')  # order of permeability from Sgem is inverse
print(perm_INV)

if isinstance(perm_INV, list):
    print("my_list is a list")
else:
    print("my_list is not a list")
    
if isinstance(perm_INV, np.ndarray):
    print("my_array is a numpy array")
else:
    print("my_array is not a numpy array")
    
def reverse(nums):
    start_index = 0
    end_index = len(nums)-1
    while end_index > start_index:
        nums[start_index],nums[end_index] = nums[end_index],nums[start_index]
        start_index = start_index + 1
        end_index = end_index -1


reverse(poro_INV)  # order of permeability from Sgem is inverse
reverse(perm_INV)  # order of permeability from Sgem is inverse

## inverting the order of the elements
#poro_INV = np.ones(len(poro_INV))*0.1
#perm_INV  = np.ones(len(poro_INV))*10

print(perm_INV)

print(max(poro_INV))
print(min(poro_INV))

class Model(DartsModel):
    def __init__(self, kvkh, n_points=128):
        # call base class constructor
        super().__init__()
        
        self.timer.node["initialization"].start()
        
        # parameters for the reservoir
        (nx, ny, nz) = (62, 62, 126)
        nb_total   = nx * ny * nz
        #Horizontal permeability was read and edited in the previous cell
        kv_kh = kvkh  # vertical horizontal permeability ratio
        perm_z = kv_kh*perm_INV
        
        # depth calculation
        depth_top = 2635 #[m]
        depth_reservoir = np.zeros(len(perm_INV))
        j=depth_top
        DZ=2.5  #[m]
        for i in range(nx*ny,nb_total+nx*ny,nx*ny):
            depth_reservoir[i-nx*ny:i] = j
            j=j+DZ
        
        # Cleaning data
        for i in range(len(perm_INV)):
            if perm_INV[i] == 0:
                perm_INV[i] = 0.1
        for i in range(len(perm_z)):
            if perm_z[i] == 0:
                perm_z[i] = 0.1
                
        for i in range(len(poro_INV)):
            if poro_INV[i] == 0:
                poro_INV[i] = 0.05
                
        
        DX= 10 # [M] all equal
        DY= 10 # [M] all equal
        self.dx = DX*np.ones(len(perm_INV)) # Units METERS
        self.dy = DY*np.ones(len(perm_INV)) # Units METERS
        self.dz= DZ*np.ones(len(perm_INV))     # Units METERS
        
        
        
        # discretize structured reservoir
        self.reservoir = StructReservoir(self.timer, nx=nx, ny=ny, nz=nz, dx=self.dx, dy=self.dy, dz=self.dz, permx=perm_INV,
                                         permy=perm_INV, permz=perm_z, poro=poro_INV, depth=depth_reservoir)
        
        # add open boundaries
        self.reservoir.set_boundary_volume(xz_minus=1e2, xz_plus=1e2, yz_minus=1e2, yz_plus=1e2,xy_minus=1e5,xy_plus=1e5) #xy_minus=1e6 to increase volume and represent overbourden
                                           
        # add well's locations
        self.jw = [30, 30]  # j well positions
        self.iw = [30, 30]  # i well positions
        # position zero is the injector and position 1 in the producer
        
        # add well
        self.reservoir.add_well("INJ")
        n_perf = 54        
        # add perforations to te payzone
        for n in range(41, n_perf):  # 2697.5 m top > 25 and 2757.5 m base > 49 then it is 50 to include 49
            self.reservoir.add_perforation(well=self.reservoir.wells[-1], i=self.iw[0], j=self.jw[0], k=n+1, 
                                       well_radius=0.16)

        # add well
        self.reservoir.add_well("PRD")
        # add perforations to te payzone        
        for n in range(70, 76):
            self.reservoir.add_perforation(self.reservoir.wells[-1], self.iw[1], self.jw[1], n+1, 0.16)
        for n in range(84, 89):
            self.reservoir.add_perforation(self.reservoir.wells[-1], self.iw[1], self.jw[1], n+1, 0.16)

        # rock heat capacity and rock thermal conduction
        hcap = np.array(self.reservoir.mesh.heat_capacity, copy=False)  #go to class
        rcond = np.array(self.reservoir.mesh.rock_cond, copy=False)
        hcap.fill(2200)
        rcond.fill(181.44)

        # create pre-defined physics for geothermal
        self.physics = Geothermal(self.timer, n_points, 1, 351, 1000, 10000, cache=False)

        # timestep parameters
        self.params.first_ts = 1e-3
        self.params.mult_ts  = 2
        self.params.max_ts   = 365

        # nonlinear and linear solver tolerance
        self.params.tolerance_newton = 1e-2
        self.params.tolerance_linear = 1e-4

        self.timer.node["initialization"].stop()
        
        
    def set_initial_conditions(self):
        # initialization with constant pressure and temperature
        #self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure=285,
        #                                            uniform_temperature=370.15)
             
        #temperature grad bar/km, Surface temperature is defined in the function
        T_grad = 28.6    # c/km
        T_surf = 10+273.15 # [K] surface temperature
        P_surf = 1 # [bar] pay attention here
        
        Top_depth = 2635 # [m]
        T_top     = T_grad*Top_depth/1000 + T_surf # estimated T reservoir to estimate initial density
        P_top     = 9.81*1000/1e5*1000*Top_depth/1000 + P_surf
        
        rho_water = 1/(_Region1(T_top, float(P_top)*0.1)['v']) #density [kg/m3], 'v' is specific volume
        p_grad = 9.81*rho_water/1e5*1000  #bar/km
        
        self.physics.set_nonuniform_initial_conditions(self.reservoir.mesh, p_grad,
                                                    T_grad,T_surf,P_surf)
        
    def set_boundary_conditions(self,T=40+273.15,Q_INJ = 0,BHP=269.4): # modified to consider injection temperature changing with time
        # activate wells with rate control for injector and producer
        for i, w in enumerate(self.reservoir.wells):
            if 'INJ' in w.name:
                w.control = self.physics.new_rate_water_inj(Q_INJ*24, T) # rate [m3/d], temperature [K], 
            else:
                w.control = self.physics.new_bhp_prod(BHP) #Controlled by bhp to guaranteed long term injection-production mass balance
                
    def export_pro_vtk(self, file_name):
        # connect to simulation array
        X = np.array(self.physics.engine.X, copy=False)
        nb = self.reservoir.mesh.n_res_blocks
        # compute temperature using pressure and enthalpy (in different units)
        temp = _Backward1_T_Ph_vec(X[0:2 * nb:2] / 10, X[1:2 * nb:2] / 18.015) # transforamtion from enthalpy to T
        # define additional arrays to the output
        local_cell_data = {'Temperature': temp,
                           'Perm': self.reservoir.global_data['permx']}
        # use export_vtk defined in the base class (DartsModel)
        self.export_vtk(file_name, local_cell_data=local_cell_data)