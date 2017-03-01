import numpy as np
import matplotlib.pyplot as plt

T100 = np.array([0.0500000, 0.1000000, 0.1500000, 0.2000000, 0.2500000, 0.3000000, 0.3500000, 0.4000000, 0.4500000, 0.5000000, 0.5500000, 0.6000000, 0.6500000, 0.7000000, 0.7500000, 0.8000000, 0.8500000, 0.9000000, 0.9500000, 1.0000000, 1.0500000, 1.1000000, 1.1500000, 1.2000000, 1.2500000, 1.3000000, 1.3500000, 1.4000000, 1.4500000, 1.5000000, 1.5500000, 1.6000000, 1.6500000, 1.7000000, 1.7500000, 1.8000000, 1.8500000, 1.9000000, 1.9500000, 2.0000000, 2.0500000, 2.1000000, 2.1500000, 2.2000000, 2.2500000, 2.3000000, 2.3500000, 2.4000000, 2.4500000, 2.5000000, 2.5500000, 2.6000000, 2.6500000, 2.7000000, 2.7500000, 2.8000000, 2.8500000, 2.9000000, 2.9500000, 3.0000000, 3.0500000, 3.1000000, 3.1500000, 3.2000000, 3.2500000, 3.3000000, 3.3500000, 3.4000000, 3.4500000, 3.5000000, 3.5500000, 3.6000000, 3.6500000, 3.7000000, 3.7500000, 3.8000000, 3.8500000, 3.9000000, 3.9500000, 4.0000000, 4.0500000, 4.1000000, 4.1500000, 4.2000000, 4.2500000, 4.3000000, 4.3500000, 4.4000000, 4.4500000, 4.5000000, 4.5500000, 4.6000000, 4.6500000, 4.7000000, 4.7500000, 4.8000000, 4.8500000, 4.9000000, 4.9500000, 5.0000000])
G100 = np.array([1.0000000, 0.9219919, 0.8600716, 0.7990608, 0.7563486, 0.7084312, 0.6741276, 0.6409509, 0.6061640, 0.5822711, 0.5620361, 0.5338788, 0.5141787, 0.4989747, 0.4841081, 0.4614621, 0.4525702, 0.4499590, 0.4269870, 0.4236950, 0.4090715, 0.3993159, 0.3973196, 0.3734209, 0.3795872, 0.3667312, 0.3628644, 0.3545760, 0.3451378, 0.3421433, 0.3328996, 0.3273625, 0.3272481, 0.3193515, 0.3151472, 0.3069817, 0.3049396, 0.3001147, 0.2937225, 0.2870814, 0.2820076, 0.2756354, 0.2757327, 0.2735447, 0.2638548, 0.2656109, 0.2573024, 0.2519913, 0.2483791, 0.2435627, 0.2422614, 0.2364898, 0.2342561, 0.2270316, 0.2239112, 0.2227129, 0.2096081, 0.2037164, 0.2032101, 0.1978246, 0.1906859, 0.1883979, 0.1817025, 0.1786994, 0.1732396, 0.1727162, 0.1616878, 0.1613589, 0.1562365, 0.1515689, 0.1494296, 0.1415673, 0.1422537, 0.1363391, 0.1319460, 0.1242153, 0.1245842, 0.1205144, 0.1088310, 0.1146941, 0.1055248, 0.1015207, 0.0997989, 0.0965013, 0.0927460, 0.0912216, 0.0877295, 0.0861136, 0.0812200, 0.0767297, 0.0754083, 0.0712012, 0.0712041, 0.0667395, 0.0648490, 0.0593148, 0.0581622, 0.0542868, 0.0563403, 0.0523934])



T250 = np.array([0.0200000, 0.0400000, 0.0600000, 0.0800000, 0.1000000, 0.1200000, 0.1400000, 0.1600000, 0.1800000, 0.2000000, 0.2200000, 0.2400000, 0.2600000, 0.2800000, 0.3000000, 0.3200000, 0.3400000, 0.3600000, 0.3800000, 0.4000000, 0.4200000, 0.4400000, 0.4600000, 0.4800000, 0.5000000, 0.5200000, 0.5400000, 0.5600000, 0.5800000, 0.6000000, 0.6200000, 0.6400000, 0.6600000, 0.6800000, 0.7000000, 0.7200000, 0.7400000, 0.7600000, 0.7800000, 0.8000000, 0.8200000, 0.8400000, 0.8600000, 0.8800000, 0.9000000, 0.9200000, 0.9400000, 0.9600000, 0.9800000, 1.0000000, 1.0200000, 1.0400000, 1.0600000, 1.0800000, 1.1000000, 1.1200000, 1.1400000, 1.1600000, 1.1800000, 1.2000000, 1.2200000, 1.2400000, 1.2600000, 1.2800000, 1.3000000, 1.3200000, 1.3400000, 1.3600000, 1.3800000, 1.4000000, 1.4200000, 1.4400000, 1.4600000, 1.4800000, 1.5000000, 1.5200000, 1.5400000, 1.5600000, 1.5800000, 1.6000000, 1.6200000, 1.6400000, 1.6600000, 1.6800000, 1.7000000, 1.7200000, 1.7400000, 1.7600000, 1.7800000, 1.8000000, 1.8200000, 1.8400000, 1.8600000, 1.8800000, 1.9000000, 1.9200000, 1.9400000, 1.9600000, 1.9800000, 2.0000000, 2.0200000, 2.0400000, 2.0600000, 2.0800000, 2.1000000, 2.1200000, 2.1400000, 2.1600000, 2.1800000, 2.2000000, 2.2200000, 2.2400000, 2.2600000, 2.2800000, 2.3000000, 2.3200000, 2.3400000, 2.3600000, 2.3800000, 2.4000000, 2.4200000, 2.4400000, 2.4600000, 2.4800000, 2.5000000, 2.5200000, 2.5400000, 2.5600000, 2.5800000, 2.6000000, 2.6200000, 2.6400000, 2.6600000, 2.6800000, 2.7000000, 2.7200000, 2.7400000, 2.7600000, 2.7800000, 2.8000000, 2.8200000, 2.8400000, 2.8600000, 2.8800000, 2.9000000, 2.9200000, 2.9400000, 2.9600000, 2.9800000, 3.0000000, 3.0200000, 3.0400000, 3.0600000, 3.0800000, 3.1000000, 3.1200000, 3.1400000, 3.1600000, 3.1800000, 3.2000000, 3.2200000, 3.2400000, 3.2600000, 3.2800000, 3.3000000, 3.3200000, 3.3400000, 3.3600000, 3.3800000, 3.4000000, 3.4200000, 3.4400000, 3.4600000, 3.4800000, 3.5000000, 3.5200000, 3.5400000, 3.5600000, 3.5800000, 3.6000000, 3.6200000, 3.6400000, 3.6600000, 3.6800000, 3.7000000, 3.7200000, 3.7400000, 3.7600000, 3.7800000, 3.8000000, 3.8200000, 3.8400000, 3.8600000, 3.8800000, 3.9000000, 3.9200000, 3.9400000, 3.9600000, 3.9800000, 4.0000000, 4.0200000, 4.0400000, 4.0600000, 4.0800000, 4.1000000, 4.1200000, 4.1400000, 4.1600000, 4.1800000, 4.2000000, 4.2200000, 4.2400000, 4.2600000, 4.2800000, 4.3000000, 4.3200000, 4.3400000, 4.3600000, 4.3800000, 4.4000000, 4.4200000, 4.4400000, 4.4600000, 4.4800000, 4.5000000, 4.5200000, 4.5400000, 4.5600000, 4.5800000, 4.6000000, 4.6200000, 4.6400000, 4.6600000, 4.6800000, 4.7000000, 4.7200000, 4.7400000, 4.7600000, 4.7800000, 4.8000000, 4.8200000, 4.8400000, 4.8600000, 4.8800000, 4.9000000, 4.9200000, 4.9400000, 4.9600000, 4.9800000, 5.0000000])
G250 = np.array([1.0000000, 0.9550518, 0.9264920, 0.9015386, 0.8719655, 0.8468234, 0.8179415, 0.8012699, 0.7680432, 0.7645700, 0.7496585, 0.7251355, 0.7033999, 0.6949266, 0.6716140, 0.6622052, 0.6492898, 0.6391980, 0.6194058, 0.6056714, 0.6041195, 0.5874617, 0.5749323, 0.5674640, 0.5527052, 0.5423358, 0.5381047, 0.5321578, 0.5205113, 0.5066076, 0.5059829, 0.5048002, 0.4903302, 0.4851996, 0.4741110, 0.4694357, 0.4634584, 0.4597631, 0.4536997, 0.4484053, 0.4429554, 0.4394323, 0.4323694, 0.4250178, 0.4199927, 0.4175634, 0.4159754, 0.4028796, 0.4054976, 0.4085682, 0.4002865, 0.3933263, 0.3828514, 0.3845394, 0.3754220, 0.3802222, 0.3775125, 0.3724847, 0.3700249, 0.3617126, 0.3668821, 0.3612823, 0.3600524, 0.3552466, 0.3552744, 0.3477201, 0.3479589, 0.3390969, 0.3388859, 0.3319063, 0.3327059, 0.3375977, 0.3327642, 0.3310123, 0.3191076, 0.3234497, 0.3237884, 0.3155317, 0.3155178, 0.3175890, 0.3173363, 0.3103623, 0.3083800, 0.3044598, 0.3086437, 0.2997429, 0.3027441, 0.3001594, 0.2946540, 0.3005175, 0.2916000, 0.2950565, 0.2922525, 0.2853395, 0.2842317, 0.2864111, 0.2780489, 0.2799840, 0.2726074, 0.2787708, 0.2731127, 0.2714497, 0.2746452, 0.2746702, 0.2590230, 0.2618937, 0.2617798, 0.2656139, 0.2591007, 0.2612912, 0.2633595, 0.2543977, 0.2517879, 0.2504775, 0.2518740, 0.2548613, 0.2534287, 0.2455718, 0.2420653, 0.2416711, 0.2445168, 0.2415989, 0.2356771, 0.2401136, 0.2360658, 0.2326426, 0.2326204, 0.2346471, 0.2363823, 0.2272455, 0.2264209, 0.2221371, 0.2246274, 0.2208072, 0.2186140, 0.2215763, 0.2188722, 0.2150270, 0.2133779, 0.2098048, 0.2056598, 0.2088414, 0.2048380, 0.2051517, 0.2039301, 0.1966562, 0.2009734, 0.1983664, 0.1959372, 0.2005264, 0.1879137, 0.1899598, 0.1870391, 0.1912452, 0.1873334, 0.1840685, 0.1884078, 0.1854400, 0.1775692, 0.1758923, 0.1762588, 0.1712725, 0.1784992, 0.1711282, 0.1656811, 0.1660892, 0.1672108, 0.1638820, 0.1622912, 0.1638515, 0.1574799, 0.1546869, 0.1543177, 0.1536292, 0.1535709, 0.1473881, 0.1520967, 0.1499895, 0.1423630, 0.1454669, 0.1423685, 0.1349225, 0.1401336, 0.1362551, 0.1393090, 0.1310773, 0.1341257, 0.1289757, 0.1303666, 0.1274959, 0.1233953, 0.1232343, 0.1216795, 0.1196889, 0.1193836, 0.1152691, 0.1172042, 0.1097082, 0.1096915, 0.1087392, 0.1119292, 0.1065737, 0.1061545, 0.1029035, 0.1042777, 0.1010683, 0.0958961, 0.0986335, 0.0963209, 0.0935140, 0.0922314, 0.0913096, 0.0930476, 0.0893801, 0.0894440, 0.0854072, 0.0858736, 0.0864678, 0.0817231, 0.0792633, 0.0804460, 0.0826476, 0.0769201, 0.0773532, 0.0766841, 0.0756125, 0.0728195, 0.0707262, 0.0696184, 0.0676806, 0.0669948, 0.0676334, 0.0660314, 0.0639520, 0.0619808, 0.0618392, 0.0586992, 0.0604483, 0.0591407, 0.0568308, 0.0584549, 0.0549290, 0.0543654, 0.0539185, 0.0524859, 0.0502760, 0.0507174, 0.0519889, 0.0493626, 0.0455812])



t = [0.01, 0.0300803, 0.0501606, 0.070241, 0.0903213, 0.1104016, 0.1304819, 0.1505622, 0.1706426, 0.1907229, 0.2108032, 0.2308835, 0.2509638, 0.2710442, 0.2911245, 0.3112048, 0.3312851, 0.3513654, 0.3714458, 0.3915261, 0.4116064, 0.4316867, 0.451767, 0.4718474, 0.4919277, 0.512008, 0.5320883, 0.5521687, 0.5722489, 0.5923293, 0.6124096, 0.63249, 0.6525703, 0.6726506, 0.6927309, 0.7128112, 0.7328916, 0.7529718, 0.7730522, 0.7931325, 0.8132128, 0.8332932, 0.8533735, 0.8734538, 0.8935341, 0.9136145, 0.9336948, 0.953775, 0.9738554, 0.9939357, 1.014016, 1.0340964, 1.0541767, 1.074257, 1.0943373, 1.1144177, 1.134498, 1.1545783, 1.1746586, 1.194739, 1.2148193, 1.2348996, 1.2549799, 1.2750602, 1.2951406, 1.3152208, 1.3353012, 1.3553815, 1.3754618, 1.3955421, 1.4156225, 1.4357028, 1.4557831, 1.4758634, 1.4959437, 1.516024, 1.5361044, 1.5561847, 1.576265, 1.5963453, 1.6164257, 1.636506, 1.6565863, 1.6766667, 1.696747, 1.7168273, 1.7369076, 1.756988, 1.7770681, 1.7971486, 1.8172289, 1.8373091, 1.8573896, 1.8774698, 1.8975501, 1.9176304, 1.9377108, 1.9577911, 1.9778714, 1.9979517, 2.0180321, 2.0381124, 2.0581927, 2.0782731, 2.0983534, 2.1184337, 2.138514, 2.1585944, 2.1786747, 2.198755, 2.2188353, 2.2389157, 2.2589958, 2.2790761, 2.2991564, 2.3192369, 2.3393172, 2.3593975, 2.3794777, 2.3995581, 2.4196384, 2.4397187, 2.4597992, 2.4798795, 2.4999597, 2.52004, 2.5401204, 2.5602007, 2.580281, 2.6003614, 2.6204418, 2.640522, 2.6606024, 2.6806827, 2.700763, 2.7208433, 2.7409237, 2.761004, 2.7810843, 2.8011647, 2.821245, 2.8413253, 2.8614056, 2.8814859, 2.9015663, 2.9216464, 2.9417267, 2.961807, 2.9818874, 3.0019677, 3.022048, 3.0421284, 3.0622087, 3.082289, 3.1023694, 3.1224497, 3.14253, 3.1626104, 3.1826907, 3.2027709, 3.2228514, 3.2429317, 3.263012, 3.2830923, 3.3031727, 3.3232529, 3.3433332, 3.3634137, 3.383494, 3.4035743, 3.4236546, 3.4437349, 3.4638152, 3.4838955, 3.5039759, 3.5240562, 3.5441364, 3.5642166, 3.5842969, 3.6043774, 3.6244577, 3.644538, 3.6646184, 3.6846987, 3.7047789, 3.7248592, 3.7449397, 3.76502, 3.7851003, 3.8051806, 3.8252609, 3.8453412, 3.8654215, 3.885502, 3.9055823, 3.9256626, 3.9457429, 3.9658232, 3.9859035, 4.0059838, 4.0260644, 4.0461445, 4.066225, 4.0863051, 4.1063856, 4.1264659, 4.1465463, 4.1666265, 4.1867069, 4.2067871, 4.2268676, 4.2469478, 4.2670283, 4.2871084, 4.307189, 4.3272691, 4.3473496, 4.3674297, 4.3875102, 4.4075904, 4.4276709, 4.447751, 4.4678315, 4.4879117, 4.5079918, 4.5280724, 4.5481524, 4.568233, 4.5883131, 4.6083937, 4.6284739, 4.6485543, 4.6686344, 4.6887149, 4.7087952, 4.7288755, 4.7489558, 4.7690363, 4.7891164, 4.8091968, 4.829277, 4.8493576, 4.8694377, 4.8895183, 4.9095984, 4.9296789, 4.949759, 4.9698395, 4.9899197, 5.0100001]
g = [0.9798003, 0.9433373, 0.9106332, 0.8808497, 0.8531358, 0.827988, 0.8042242, 0.7819953, 0.7611346, 0.7417858, 0.7236349, 0.7063062, 0.6894959, 0.673915, 0.6592969, 0.6452666, 0.631905, 0.6189307, 0.6075261, 0.5956507, 0.5841034, 0.5735535, 0.5641173, 0.5534044, 0.5444099, 0.5354072, 0.5285558, 0.5182548, 0.5109252, 0.5046207, 0.4964894, 0.488297, 0.4828015, 0.4756241, 0.4674257, 0.4635856, 0.4561774, 0.4514218, 0.444792, 0.4405735, 0.434819, 0.4279687, 0.4259177, 0.4180259, 0.4155963, 0.4092824, 0.4064189, 0.4017211, 0.3992741, 0.3956648, 0.3908532, 0.3848906, 0.3816059, 0.377251, 0.3737123, 0.3732704, 0.3670268, 0.3652869, 0.3615445, 0.3598042, 0.3563428, 0.3537503, 0.3520228, 0.3455104, 0.3449193, 0.3401729, 0.339383, 0.3370182, 0.3340539, 0.3309914, 0.3302265, 0.3291696, 0.3239483, 0.3214822, 0.3200881, 0.3175776, 0.3205604, 0.3140031, 0.3114183, 0.3112214, 0.3042363, 0.3055045, 0.3037177, 0.301093, 0.2964551, 0.2990299, 0.2909793, 0.2911953, 0.2895733, 0.2857277, 0.2850147, 0.2890508, 0.286876, 0.2825479, 0.2790656, 0.2831613, 0.2799424, 0.2775154, 0.272472, 0.2698035, 0.2661871, 0.2694947, 0.2684036, 0.2649558, 0.2585015, 0.2607577, 0.2607856, 0.2597452, 0.2590155, 0.2513425, 0.2557385, 0.25539, 0.247287, 0.2497781, 0.2471318, 0.2441113, 0.2417419, 0.2451208, 0.2388599, 0.235296, 0.2340182, 0.2363589, 0.2346767, 0.2359747, 0.226886, 0.2280738, 0.2263902, 0.224368, 0.2269911, 0.2207875, 0.2184601, 0.2189817, 0.2190477, 0.2145169, 0.2121933, 0.2088551, 0.2121766, 0.2087067, 0.2069482, 0.2055303, 0.2059401, 0.208197, 0.1993071, 0.201458, 0.1977797, 0.191149, 0.1933481, 0.192562, 0.1928081, 0.1915811, 0.1878316, 0.1900191, 0.1806498, 0.1827671, 0.1788121, 0.1773146, 0.1830969, 0.1765616, 0.1718066, 0.1676733, 0.1705097, 0.1652724, 0.1631669, 0.1571055, 0.1643511, 0.1616818, 0.1580416, 0.152401, 0.163261, 0.1545914, 0.1550673, 0.1514083, 0.1417249, 0.148499, 0.1422588, 0.1434121, 0.1366974, 0.1403497, 0.1325942, 0.1364492, 0.1380059, 0.1397359, 0.126273, 0.1254442, 0.1260027, 0.1239855, 0.1240504, 0.1191578, 0.1205528, 0.1098035, 0.1191182, 0.1205616, 0.1211653, 0.1159217, 0.1081335, 0.1175955, 0.1090865, 0.108327, 0.1063523, 0.0987173, 0.1026108, 0.1032383, 0.104086, 0.1038375, 0.0943245, 0.0981819, 0.0917449, 0.095769, 0.0938004, 0.0860774, 0.0852963, 0.0927917, 0.0825183, 0.0789799, 0.0844813, 0.0830605, 0.0815756, 0.0732507, 0.0765095, 0.0809513, 0.080282, 0.0781006, 0.0715543, 0.0741776, 0.0716124, 0.0692019, 0.0691616, 0.0681633, 0.0665342, 0.0704983, 0.0701031, 0.0711707, 0.0626267, 0.0597334, 0.0617144, 0.060469, 0.0614918, 0.0570772, 0.0600805, 0.0559629, 0.0531729, 0.0524054, 0.0521472, 0.0524737, 0.0493251, 0.0462201, 0.0535328, 0.0480334, 0.0493713, 0.0426921]


plt.plot(t, g, 'r')
# plt.plot(T100 - 0.025, G100*0.952, 'b')
plt.plot(T250 - T250[0]/2, G250*0.980, 'g')
plt.show()