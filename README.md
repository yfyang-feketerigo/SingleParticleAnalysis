# SingleParticleAnalysis

用于计算KA模型在剪切过程中的单粒子动力学、静态结构等分析。

原始数据基于LAMMPS输出。目前使用data文件(LAMMPS write_data命令产生)作为输入文件。

-----
PaAn_Settings.json为配置文件：

  `"_comment": "路径名记得添加 '/'",`
  
  `"wi": 10,`                                
  剪切速率，使用无量纲魏森贝格数  
  
  `"tau_alpha": 12000,`                       
  α松弛时间，LJ单位(原则上可以使用其他单位制，应与LAMMPS使用单位保持一致）  
  
  `"data_fpath": "./",`  
  data文件存放路径
  
  `"fname_prefix": "data.restart.shear.wi.10.",`     
  data文件前缀，后缀应为时间步，例如data.restart.shear.wi.10.240000
  
  `"shear_data_pairstyle": "single",`  
  data文件包含势能信息的格式与形式，目前包括
  
  single: 存储相同粒子间的势能信息，不同粒子间势能信息丢失(LAMMPS可能会按照特定算法恢复，但是往往并非正确结果，需要注意)
  
  pair: 存储相同以及不同粒子间的势能信息
  
  none：没有势能信息包含在data文件中
  
  `"fname_postfix": null,`                           
  文件后缀，用.分割，应对时间步在中间的情况。备用方案，未经严格测试。
  
  `"delta_step": 240000,`                            
  处理间隔，即每隔delta_step步进行一次处理计算。单位为MD时间步，并非时间单位。
  
  `"moment_number": 1,`  
  总时刻数，即总共需要计算多少个瞬间。
  
  `"start_moment": 1,`                               
  起始时刻编号，start_moment * delta_step 即为第一次计算的时间步
  
  `"equi_data_fpath": "./",`                         
  平衡态构象文件存放路径
  
  `"equi_config_fname": "data.T0.42.Equi",`          
  平衡态构象文件名称
  
  `"equi_data_pairstyle": "single",`                 
  平衡态构象文件包含的势能信息与形式。
  
  `"output_path": "./test/",`                        
  输出文件路径
  
  `"CN_rcut": 1.0,`                                  
  配位数计算的截断半径
  
  `"computeCN": true,`                               
  计算配位数
  
  `"computeMSD": true,`                              
  计算MSD
  
  `"computeMSDnonAffine": true,`                     
  计算非仿射MSD
  
  `"MSDnonAffine_ave_gradient": true,`               
  计算梯度方向时间平均后的MSD(物理意义上被废弃的功能，但是没删)
  
  `"MSDnonAffine_t0": true,`                         
  计算按t0时刻作为梯度位置的MSD(物理意义上被废弃的功能，但是没删)
  
  `"computeFlowDisplacement": true,`  
  计算在delta_step时长内，每个粒子在流场方向上的位移
  
  `"computeFlowAveVelocity": true`                  
  计算在delta_step时长内，每个粒子在流场方向上的平均速度  

---
包括如下文件：  
`configuration.h`,`configuration.cpp`：构象类，将data文件作为构象读入  

`input.h`,`input.cpp`：Input类，按行读入数据文件

`particle.h`,`particle.cpp`: Particle结构

`main.cpp`:主程序


  
