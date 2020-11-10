from __future__ import absolute_import, division, print_function

import sys, math
from cctbx.array_family import flex
from six.moves import range
from six.moves import zip
from libtbx.utils import null_out
from iotbx.map_model_manager import map_model_manager

# Binned rms amplitude for beta-galactosidase calculated from 1 A - inf
# with k_sol=0  (no solvent)

rms_fc_list = """
d_min   rmsFc
14.8935 27428.7
8.7965  9656.1
7.3874  8589.7
6.5955  8439.5
6.0632  8732.1
5.6666  9441.4
5.3566  9666.1
5.1085  10218.7
4.897   10634.3
4.7177  10939.9
4.5664  10926.3
4.4289  10606.2
4.3077  10309.6
4.1972  9365.6
4.1     8950.7
4.0098  8536.7
3.9278  8025.5
3.8517  8133.3
3.7793  7798.6
3.7141  7398
3.6529  7173.5
3.5949  6839.8
3.5413  6978
3.4908  6838
3.4416  6525.5
3.3965  6349.7
3.3533  6226.1
3.3119  6216
3.2726  6040.4
3.2348  6097.3
3.1994  5721.7
3.1652  5771.8
3.1317  5756.4
3.101   5480.2
3.0702  5583.7
3.0421  5317.7
3.0136  5460.2
2.9864  5207.4
2.9606  5259.4
2.9356  5089.8
2.9111  5273.6
2.8873  5087
2.8645  5111.6
2.8425  5260.2
2.8208  4967.6
2.8002  5052.1
2.7796  5086.5
2.7601  5217.9
2.741   4970.3
2.7226  5013.3
2.7045  5067.9
2.6865  5120.5
2.6693  4810.8
2.6529  5161.3
2.6364  5038.5
2.6206  4812.4
2.6046  4980.8
2.5899  4894.3
2.5752  4907.6
2.5607  5007.3
2.5465  5057.8
2.5323  5053.6
2.5185  4956.8
2.5056  5037.1
2.4928  4719.8
2.4798  5245.1
2.4672  4924.7
2.455   4913.6
2.4431  4872.6
2.4313  4843.4
2.4197  4934.3
2.4085  4828.5
2.3974  4837.6
2.3864  5012.5
2.3756  4823.7
2.3652  4851
2.3548  4848.3
2.3445  4925
2.3344  4794.6
2.3249  5099.7
2.3149  4886.7
2.3058  4938.6
2.2963  4831.9
2.2869  4932.6
2.2778  4955
2.269   4834
2.2603  4922.9
2.2515  4629
2.243   4956.5
2.2346  4622.5
2.2263  4760.7
2.2185  4681.7
2.2103  4740.7
2.2024  4687.5
2.1945  4845.7
2.1868  4684.2
2.1793  4711.9
2.1719  4687.6
2.1646  4789.3
2.1571  4857.2
2.1499  4775
2.1429  4578.3
2.136   4638.5
2.1291  4589.6
2.1222  4639.6
2.1154  4688.8
2.1088  4538.1
2.1023  4759.4
2.0957  4541.6
2.0893  4467.1
2.0829  4494.1
2.0769  4556.1
2.0707  4522.1
2.0646  4391.1
2.0587  4372.9
2.0526  4492.9
2.0466  4369.9
2.0408  4305.6
2.0351  4429.5
2.0294  4311.8
2.0239  4458.9
2.0182  4340.4
2.0126  4479.6
2.0073  4293
2.0019  4444.5
1.9966  4024.8
1.9913  4402.5
1.9861  4320.9
1.9809  4264.3
1.9758  4272.1
1.9707  4193.8
1.9657  4306.6
1.9609  4192.2
1.9559  4224
1.951   4105.9
1.9462  3959.2
1.9413  4212.7
1.9366  4207.6
1.9321  4119.3
1.9274  3971.2
1.9227  3997.4
1.9183  4088
1.9137  4017.2
1.9093  4027.2
1.9048  3956.5
1.9004  3789.9
1.8962  3892.8
1.8919  3960
1.8877  3801.2
1.8834  3710.2
1.8791  3909.5
1.8752  3969.3
1.871   3778.8
1.867   3806.7
1.8628  3664.6
1.8589  3777.6
1.8549  3754.3
1.8509  3740
1.8471  3669.2
1.8432  3591.1
1.8394  3707.7
1.8355  3570.3
1.8318  3626.8
1.8281  3601.8
1.8243  3585.9
1.8207  3595.2
1.817   3704.4
1.8134  3664.5
1.8098  3442.9
1.8062  3517.1
1.8027  3507.4
1.7992  3522.2
1.7956  3543.2
1.7923  3477.4
1.7889  3553.9
1.7854  3459.4
1.782   3488.1
1.7787  3489.7
1.7753  3329.4
1.7721  3591.5
1.7688  3320.3
1.7656  3352.5
1.7623  3365.4
1.7591  3370.3
1.7559  3425.1
1.7527  3332.2
1.7495  3309.9
1.7465  3350.8
1.7433  3408
1.7403  3342.2
1.7372  3306.2
1.7343  3275.2
1.7312  3367.8
1.7282  3301.3
1.7252  3216.8
1.7223  3385.3
1.7193  3231.7
1.7165  3230.7
1.7137  3379.6
1.7108  3270.5
1.7078  3375.6
1.705   3328.2
1.7022  3040.2
1.6994  3256.9
1.6966  3125.6
1.6939  3225.9
1.6911  3186.4
1.6884  3206.2
1.6857  3163.7
1.683   3151.8
1.6804  3162.4
1.6777  3150.6
1.6751  3103.6
1.6725  3096.7
1.6699  3117.2
1.6673  3136.9
1.6647  3081.4
1.6622  3044.6
1.6596  3034.3
1.657   3017.7
1.6546  2990.9
1.6521  3068.2
1.6496  3066.4
1.6471  2975.2
1.6447  2970.8
1.6423  3040.8
1.6399  3072.4
1.6374  2919.9
1.6351  3023.2
1.6327  3071.4
1.6303  3005.4
1.6279  2933.8
1.6256  3095.9
1.6233  3045.8
1.621   2909.1
1.6187  2930.7
1.6164  3036.3
1.6141  2937.7
1.6119  2938.1
1.6097  2953.6
1.6074  3058.2
1.6051  2904.5
1.6029  2953.9
1.6008  2823.3
1.5986  3006.7
1.5965  2942.5
1.5943  3016
1.5921  2916.4
1.59    3037.2
1.5878  3006.2
1.5857  2819.5
1.5836  2896.3
1.5815  2828.3
1.5795  2819.2
1.5773  2872.6
1.5753  3014.2
1.5733  2907.1
1.5712  2809
1.5691  2910
1.5672  2752.1
1.5651  3028.9
1.5631  2778.4
1.5612  2803.5
1.5592  2845.2
1.5573  2919.3
1.5553  2885
1.5533  2805.4
1.5514  2833.2
1.5494  2851.3
1.5475  2737.8
1.5456  2865.3
1.5438  2876.5
1.5418  2815.8
1.54    2717.8
1.5381  2785.3
1.5362  2786.2
1.5344  2911.6
1.5325  2770.4
1.5307  2848.2
1.5289  2785.3
1.527   2760.1
1.5252  2698.8
1.5234  2886.7
1.5216  2820.9
1.5199  2703.3
1.5181  2721.4
1.5163  2663.4
1.5145  2761.1
1.5128  2787.6
1.511   2814.4
1.5093  2821.9
1.5076  2732.4
1.5059  2715.3
1.5042  2806.9
1.5025  2734.3
1.5008  2849.6
1.4991  2671.9
1.4974  2671.6
1.4957  2830.2
1.4941  2703.5
1.4924  2797.7
1.4907  2711.3
1.4891  2751.4
1.4875  2850.9
1.4858  2750.8
1.4842  2575.8
1.4826  2731.3
1.4809  2837.1
1.4793  2763.9
1.4778  2697.9
1.4761  2785.9
1.4746  2686.2
1.473   2719
1.4714  2672.9
1.4699  2829.7
1.4683  2807.6
1.4668  2685.9
1.4653  2672.9
1.4637  2642.8
1.4622  2815
1.4606  2531.2
1.4592  2728.3
1.4577  2725.9
1.4561  2752.1
1.4547  2723.1
1.4531  2710
1.4517  2713.8
1.4502  2712.2
1.4487  2589.5
1.4472  2636.1
1.4457  2755.8
1.4444  2637.9
1.4429  2606
1.4414  2704.1
1.44    2690.4
1.4386  2722.4
1.4371  2599.4
1.4358  2731.4
1.4343  2628.3
1.4329  2707.3
1.4315  2662.9
1.4301  2769.2
1.4287  2605.1
1.4273  2542
1.4259  2689.2
1.4245  2608.5
1.4232  2638.8
1.4218  2630.3
1.4205  2686.3
1.4191  2728.7
1.4178  2589.5
1.4164  2681.2
1.415   2650.8
1.4137  2606.4
1.4124  2711.1
1.4111  2729.2
1.4097  2563.5
1.4084  2671.9
1.4071  2645.3
1.4058  2642.2
1.4045  2745.9
1.4032  2679
1.4019  2746.6
1.4006  2642.5
1.3993  2634.9
1.3981  2615.7
1.3968  2686
1.3955  2604.9
1.3943  2688.5
1.393   2708.4
1.3918  2647.7
1.3905  2635.1
1.3892  2584.1
1.388   2524.6
1.3868  2767.5
1.3856  2769.5
1.3843  2528.1
1.3831  2709.6
1.3818  2670.2
1.3807  2533.2
1.3795  2599.6
1.3782  2587.2
1.377   2654.2
1.3758  2619.5
1.3746  2608.1
1.3735  2618.2
1.3723  2669.2
1.3711  2731.1
1.3699  2449.7
1.3687  2624.4
1.3676  2667.5
1.3664  2577.7
1.3653  2548.4
1.3641  2687.9
1.3629  2638.9
1.3618  2616.6
1.3607  2545.5
1.3595  2515.1
1.3584  2552.9
1.3573  2484.9
1.3561  2638.1
1.355   2656.9
1.3539  2468.9
1.3527  2563.8
1.3516  2685.5
1.3505  2622.2
1.3494  2591.8
1.3483  2645.3
1.3472  2580.5
1.3461  2525.2
1.345   2549.1
1.3439  2530.4
1.3428  2648.4
1.3417  2545.7
1.3407  2644
1.3396  2646.4
1.3385  2511.4
1.3374  2670.6
1.3364  2575.9
1.3353  2573.5
1.3343  2526.2
1.3332  2585
1.3322  2531
1.3311  2589.2
1.33    2563.7
1.329   2637.5
1.328   2499.9
1.3269  2551.8
1.3259  2562.9
1.3249  2641.1
1.3238  2636.9
1.3228  2606.8
1.3218  2635
1.3208  2497.1
1.3198  2479.9
1.3188  2625.3
1.3178  2602.3
1.3168  2608.4
1.3157  2668.2
1.3147  2575.8
1.3137  2555
1.3128  2542.9
1.3118  2521.9
1.3108  2592.6
1.3098  2574.2
1.3088  2668.4
1.3078  2552.5
1.3069  2588.9
1.3059  2578.1
1.3049  2560.1
1.3039  2565.2
1.303   2571
1.302   2509.6
1.301   2515.5
1.3001  2613.3
1.2992  2593.9
1.2982  2495.1
1.2973  2442.2
1.2963  2587.1
1.2954  2629.3
1.2945  2494.6
1.2935  2588.5
1.2926  2563.5
1.2916  2493.9
1.2907  2685.4
1.2898  2585.1
1.2889  2526.5
1.288   2557.2
1.287   2485.5
1.2861  2621.4
1.2852  2561.2
1.2843  2505
1.2834  2627.7
1.2825  2476.7
1.2816  2558.4
1.2807  2537.7
1.2798  2591.6
1.2789  2618.9
1.278   2500.3
1.2771  2537.7
1.2762  2556.8
1.2754  2545.5
1.2745  2604.9
1.2736  2568
1.2727  2565.1
1.2718  2393.1
1.2709  2639.1
1.2701  2639
1.2692  2504.1
1.2684  2567.8
1.2675  2546.4
1.2666  2548.4
1.2658  2559
1.2649  2571.1
1.2641  2546.4
1.2632  2528.5
1.2624  2526.8
1.2615  2553.9
1.2607  2475.2
1.2598  2561.4
1.259   2565.5
1.2581  2509.1
1.2573  2511.7
1.2565  2609.8
1.2556  2482.8
1.2548  2560.8
1.254   2490.3
1.2532  2482.6
1.2524  2586.7
1.2515  2529.6
1.2507  2475
1.2499  2623.4
1.2491  2541.3
1.2483  2503.7
1.2475  2560.6
1.2467  2524.9
1.2459  2575.7
1.2451  2516.9
1.2443  2463.2
1.2434  2504.6
1.2427  2654.8
1.2419  2519
1.2411  2470.4
1.2403  2518.7
1.2395  2584.5
1.2387  2614.3
1.2379  2571.4
1.2372  2534.8
1.2364  2567.5
1.2356  2440.5
1.2348  2448.3
1.234   2659
1.2333  2437.7
1.2325  2504.7
1.2317  2560
1.2309  2621.9
1.2302  2558.5
1.2294  2400.1
1.2286  2572.5
1.2279  2507.4
1.2272  2430.2
1.2264  2539.4
1.2256  2555.1
1.2249  2543.6
1.2241  2524.4
1.2234  2677.1
1.2226  2538.8
1.2219  2472.2
1.2211  2617.8
1.2204  2445.6
1.2197  2494.7
1.2189  2593.7
1.2182  2514.1
1.2175  2517.6
1.2167  2579.7
1.216   2587.4
1.2153  2427.5
1.2145  2545.1
1.2138  2499.2
1.2131  2624.3
1.2124  2453.6
1.2116  2552
1.2109  2557.1
1.2102  2516.5
1.2095  2495.7
1.2088  2520.6
1.2081  2550.2
1.2074  2467.5
1.2067  2480.4
1.206   2581
1.2052  2571.5
1.2045  2460.8
1.2039  2560.4
1.2031  2505.6
1.2024  2563.7
1.2017  2579.5
1.201   2529.6
1.2004  2475.8
1.1997  2558.8
1.199   2530.2
1.1983  2553.2
1.1976  2519.2
1.1969  2499.2
1.1962  2583.7
1.1955  2492.7
1.1949  2590.1
1.1942  2476.8
1.1935  2458.8
1.1928  2527.4
1.1922  2536.9
1.1915  2576.7
1.1908  2437.7
1.1901  2480.7
1.1895  2601
1.1888  2528.4
1.1881  2622.2
1.1875  2519.4
1.1868  2546.6
1.1861  2613.2
1.1855  2516.3
1.1848  2371.9
1.1842  2571.1
1.1835  2525.4
1.1829  2481.9
1.1822  2507.7
1.1816  2510.8
1.1809  2542.3
1.1803  2560.3
1.1796  2528
1.179   2537
1.1783  2411.3
1.1777  2520.2
1.177   2554.6
1.1764  2468.4
1.1758  2485.9
1.1751  2584.6
1.1745  2487.6
1.1738  2600
1.1732  2522.2
1.1726  2530.8
1.172   2669
1.1713  2473.5
1.1707  2507.8
1.1701  2492.3
1.1694  2356.6
1.1688  2594.2
1.1682  2558.9
1.1676  2387.9
1.167   2558.6
1.1663  2524.4
1.1657  2556.6
1.1651  2611.1
1.1645  2509.1
1.1639  2457.6
1.1633  2578
1.1627  2650.5
1.1621  2407.5
1.1615  2507.3
1.1608  2404.1
1.1602  2540.5
1.1596  2467.9
1.159   2418.2
1.1584  2484.1
1.1578  2558.3
1.1572  2432.4
1.1566  2555
1.156   2512.2
1.1554  2465.6
1.1548  2567.6
1.1542  2512.8
1.1537  2551.8
1.1531  2527.5
1.1525  2576.5
1.1519  2394.6
1.1513  2468.9
1.1507  2508.9
1.1501  2624.5
1.1495  2597.2
1.1489  2392.1
1.1484  2516.8
1.1478  2526.2
1.1472  2453.9
1.1467  2562.5
1.1461  2543.5
1.1455  2490.4
1.1449  2480.3
1.1443  2482.4
1.1438  2492.9
1.1432  2513.9
1.1426  2519.7
1.1421  2517.1
1.1415  2422.5
1.1409  2440.6
1.1404  2560.7
1.1398  2584
1.1392  2573.1
1.1387  2454
1.1381  2453.5
1.1376  2591.1
1.137   2490
1.1364  2467.5
1.1359  2454.5
1.1353  2578.7
1.1348  2528.4
1.1342  2422.8
1.1337  2562.1
1.1331  2480.5
1.1326  2544.9
1.132   2496.6
1.1315  2482.9
1.1309  2458.8
1.1304  2457
1.1298  2444.8
1.1293  2386.7
1.1288  2536.1
1.1282  2529.2
1.1277  2477.7
1.1271  2461.6
1.1266  2494.9
1.1261  2507.8
1.1255  2478
1.125   2451.5
1.1245  2436.3
1.1239  2368.3
1.1234  2424.9
1.1229  2431.1
1.1223  2521.9
1.1218  2470.9
1.1213  2460.5
1.1207  2516.8
1.1202  2475.5
1.1197  2405.4
1.1192  2439
1.1186  2441
1.1181  2496.2
1.1176  2446.5
1.1171  2330.1
1.1166  2475.1
1.116   2415.3
1.1155  2480.8
1.115   2350
1.1145  2333.8
1.114   2451.5
1.1135  2326.5
1.113   2340.7
1.1124  2489.9
1.1119  2588.4
1.1114  2392.3
1.1109  2361
1.1104  2463.2
1.1099  2437.4
1.1094  2378.8
1.1089  2395.3
1.1084  2387.1
1.1079  2367.3
1.1074  2405.4
1.1069  2374.7
1.1064  2255.9
1.1059  2446.1
1.1054  2405.4
1.1049  2516.1
1.1044  2350.1
1.1039  2480.7
1.1034  2356.2
1.1029  2326.6
1.1024  2395.5
1.1019  2344.2
1.1014  2496.7
1.1009  2332.3
1.1004  2494.9
1.1     2320
1.0995  2390.3
1.099   2377.8
1.0985  2368.8
1.098   2425.6
1.0975  2404.4
1.097   2336.8
1.0966  2464.6
1.0961  2356.6
1.0956  2430.1
1.0951  2284.9
1.0946  2427
1.0942  2386.6
1.0937  2374.1
1.0932  2294.3
1.0927  2355.8
1.0923  2439.6
1.0918  2422.4
1.0913  2226.2
1.0908  2281.2
1.0903  2320.4
1.0899  2427.3
1.0894  2357.9
1.0889  2417.8
1.0885  2316.3
1.088   2363.3
1.0875  2290.5
1.0871  2464.4
1.0866  2299.1
1.0861  2374.1
1.0857  2304.2
1.0852  2315.3
1.0848  2351.3
1.0843  2300.2
1.0838  2324.1
1.0834  2334
1.0829  2263
1.0825  2315.3
1.082   2334.8
1.0815  2345.9
1.0811  2264.1
1.0806  2274
1.0802  2403.2
1.0797  2287
1.0793  2312.9
1.0788  2178.6
1.0784  2341.7
1.0779  2301.3
1.0775  2251.4
1.077   2234.9
1.0766  2273.2
1.0761  2356.6
1.0757  2277.6
1.0752  2215
1.0748  2230.5
1.0743  2271.5
1.0739  2204.3
1.0734  2284.8
1.073   2287.5
1.0726  2259.1
1.0721  2259.2
1.0717  2207.4
1.0712  2213.1
1.0708  2232.9
1.0704  2187.6
1.0699  2312.8
1.0695  2414.7
1.069   2190.7
1.0686  2210.9
1.0682  2188.3
1.0677  2219.1
1.0673  2218.6
1.0669  2181.9
1.0664  2144.1
1.066   2286
1.0656  2194.8
1.0652  2248.8
1.0647  2174.8
1.0643  2185.5
1.0639  2096.8
1.0634  2247.4
1.063   2249.1
1.0626  2164.5
1.0622  2162.2
1.0617  2188.6
1.0613  2181.3
1.0609  2204.4
1.0605  2183.7
1.06    2202.1
1.0596  2203.2
1.0592  2176
1.0588  2172.2
1.0584  2212.9
1.0579  2239.1
1.0575  2057.6
1.0571  2170
1.0567  2146.5
1.0563  2134.2
1.0559  2181.7
1.0554  2181.9
1.055   2180.9
1.0546  2160.3
1.0542  2124.1
1.0538  2115.9
1.0534  2158.2
1.053   2068.1
1.0526  2118.6
1.0522  2140.4
1.0518  2083.1
1.0513  2108
1.0509  2137.2
1.0505  2123.3
1.0501  2115.2
1.0497  2164.7
1.0493  2038.6
1.0489  2134.2
1.0485  2019.4
1.0481  2174.9
1.0477  2076.7
1.0473  2134.1
1.0469  2034.2
1.0465  2094.4
1.0461  2082.2
1.0457  2120.5
1.0453  2015.8
1.0449  2089.4
1.0445  2073.7
1.0441  2115.3
1.0437  1993.2
1.0433  2117.1
1.0429  2067
1.0425  2051.2
1.0421  2008.4
1.0417  2036
1.0413  2118.7
1.041   2031.3
1.0406  2081.8
1.0402  2083.1
1.0398  2084.4
1.0394  1964.5
1.039   2013.8
1.0386  2010.1
1.0382  1971.5
1.0378  2044
1.0374  1964.5
1.0371  1985.2
1.0367  2044.1
1.0363  2065.4
1.0359  1980.4
1.0355  2026.5
1.0351  1993.5
1.0348  2062.8
1.0344  1927
1.034   1956.2
1.0336  1910.5
1.0332  1948.3
1.0329  2017.1
1.0325  2020.7
1.0321  1916.6
1.0317  1926.7
1.0313  1920.9
1.031   2015.1
1.0306  1923.3
1.0302  1962.7
1.0298  2004.8
1.0294  1967.4
1.0291  1987.2
1.0287  1978.1
1.0283  1907.8
1.028   1934.1
1.0276  1893.9
1.0272  1910.7
1.0268  1940.4
1.0265  1889.2
1.0261  1920
1.0257  1986.7
1.0254  1951.8
1.025   1901.3
1.0246  1910.8
1.0243  1928.7
1.0239  1849.4
1.0235  1893.2
1.0232  1871.1
1.0228  1934
1.0224  1926.5
1.0221  1782.2
1.0217  1895.1
1.0213  1838
1.021   1840.8
1.0206  1814.7
1.0202  1875.3
1.0199  1919.6
1.0195  1804.4
1.0192  1807.9
1.0188  1865.9
1.0184  1839.4
1.0181  1835
1.0177  1823.3
1.0174  1776.3
1.017   1767
1.0167  1847.1
1.0163  1838.4
1.0159  1738.2
1.0156  1785.1
1.0152  1809
1.0149  1689.4
1.0145  1839.2
1.0142  1795.6
1.0138  1704
1.0135  1762.3
1.0131  1726.4
1.0128  1751.6
1.0124  1702.1
1.0121  1706.2
1.0117  1748.2
1.0114  1695.8
1.011   1678
1.0107  1719.9
1.0103  1580.7
1.01    1746.7
1.0096  1685.2
1.0093  1637.4
1.0089  1594.6
1.0086  1592.8
1.0083  1581.1
1.0079  1491.4
1.0076  1526.5
1.0072  1506.9
1.0069  1528.4
1.0065  1507.8
1.0062  1504.3
1.0058  1438.5
1.0055  1458.2
1.0052  1411.7
1.0048  1424.4
1.0045  1437.9
1.0041  1365.2
1.0038  1373.5
1.0035  1346.5
1.0031  1341.3
1.0028  1294.2
1.0024  1334.1
1.0021  1351
1.0018  1291.3
1.0014  1333.7
1.0011  1196.1
1.0008  1269.1
1.0004  1245.6
1.0001  1213.5
0.9998  1206.6
"""

def get_expected_ssqr_list(
      n_bins = None,
      d_min = None, # minimum resolution of data to use
      expected_ssqr_list_rms = None,
      map_model_manager = None,
      out = sys.stdout):

  assert n_bins is not None
  assert d_min is not None
  assert expected_ssqr_list_rms is not None

  if str(expected_ssqr_list_rms) == 'Auto':
    expected_ssqr_list_rms = 0.25 * d_min

  # Estimate E**2 vs resolution based on rms error expected_ssqr_list_rms
  # Assume constant total error, amplitudes fall off exp(-b_eff sthol2)
  #  b_eff = (8 * 3.14159**2/3.) * expected_ssqr_list_rms**2
  #  sigmaa = exp(-b_eff sthol2)
  #  E**2 = 1/(1-sigmaa**2)

  map_coeffs = map_model_manager.get_any_map_manager(
     ).map_as_fourier_coefficients(d_min = d_min)
  from iotbx.map_model_manager import get_map_coeffs_as_fp_phi
  f_array_info = get_map_coeffs_as_fp_phi(
    map_coeffs, d_min= d_min, n_bins = n_bins)
  f_array = f_array_info.f_array
  # Approx fall-off with resolution based on rms in A ...as in sigmaa
  b_eff = (8 * 3.14159**2/3.) * expected_ssqr_list_rms**2

  dsd = f_array.d_spacings().data()
  ssqr_list = flex.double()
  for i_bin in f_array.binner().range_used():
    sel       = f_array.binner().selection(i_bin)
    d         = dsd.select(sel)
    d_avg     = flex.mean(d)

    sthol2 = 0.25/d_avg**2
    sigmaa = math.exp(max(-20.,min(20., -b_eff * sthol2)))
    ssqr = (1-sigmaa**2)/max(1.e-10,sigmaa**2)
    ssqr_list.append(ssqr)
  return ssqr_list

class approx_amplitude_vs_resolution:
  def __init__(self,
      n_bins = 1000,
      d_min = None, # minimum resolution of data to use
      resolution = None,  # nominal resolution
      map_model_manager = None,
      model = None,
      generate_mock_rms_fc_list = True,
      k_sol = None,
      b_sol = None,
      find_k_sol_b_sol = True,
      map_id_for_optimization = 'map_manager',
      out = sys.stdout):
    #  If d_min and map_model_manager supplied and generate_mock_rms_fc_list,
    #  generate a  model-based map in that cell and calculate binned values
    # If model, use it to generate fc instead of mock.  Set B-iso = 0

    self.rms_fc_list = None
    if generate_mock_rms_fc_list:
      if model:
        model = model.deep_copy()
        model.set_b_iso(flex.double(model.get_b_iso().size(),0.)) # XXX
        u_cart=model.get_xray_structure().scatterers().extract_u_cart(model.get_xray_structure().crystal_symmetry().unit_cell())
        model.get_xray_structure().set_u_cart(u_cart)
        model.set_xray_structure(model.get_xray_structure())

      assert d_min and map_model_manager
      if find_k_sol_b_sol and map_model_manager.get_map_manager_by_id(
         map_id_for_optimization):
        from cctbx.maptbx.refine_sharpening import get_effective_b
        rms_fo_info = map_model_manager.get_rms_f_list(
          map_id = map_id_for_optimization,
          n_bins = n_bins,
          resolution = resolution,
          d_min = d_min)
        n_bins_use = rms_fo_info.n_bins_use

        # Figure out k_sol and b_sol so that rms_f vs resolution from model
        #  is related by simple Wilson B to rms_f vs resolution from map
        #  Use grid search to not take very long and get approx answer
        if (k_sol is None) and (b_sol is None):
          kb_list= [ [0,0], [0.1,20], [0.1,50],
                      [0.2,20], [0.2,50],
                      [0.3,20], [0.3,50],
                   ]

          best_rms = None
          best_kb = None
          best_rms_fc_list = None
          for k_sol,b_sol in kb_list:
            self.generate_mock_rms_fc_list(n_bins, d_min , map_model_manager,
              model = model,
              k_sol = k_sol,
              b_sol = b_sol,)

            ratio_list = rms_fo_info.rms_f_list/self.rms_fc_list
            info = get_effective_b(values = ratio_list,
              sthol2_values = rms_fo_info.sthol2_list)
            scaled_fc_list = self.rms_fc_list * info.calc_values
            rms = ((flex.pow2(rms_fo_info.rms_f_list[:n_bins_use] -
              scaled_fc_list[:n_bins_use])).min_max_mean().mean)**0.5
            if best_rms is None or rms < best_rms:
              best_rms_fc_list = rms_fc_list
              best_rms = rms
              best_kb = [k_sol,b_sol]
          k_sol, b_sol = best_kb
          if len(kb_list) > 1:  # redo calculation
            self.generate_mock_rms_fc_list(n_bins, d_min , map_model_manager,
              model = model,
              k_sol = k_sol,
              b_sol = b_sol,)
          print("Overall approximate B-value: "+
             "%.2f A**2 \n(rms = %.1f, k_sol=%.2f  b_sol=%.2f)" %(
             info.effective_b, rms,k_sol,b_sol,
              ), file = out)
          self.k_sol = k_sol
          self.b_sol = b_sol
          self.rms_fc_list = best_rms_fc_list

      # Othewise just generate the list
      self.generate_mock_rms_fc_list(n_bins, d_min , map_model_manager,
        model = model,
        k_sol = k_sol,
        b_sol = b_sol,
        n_bins_use = n_bins_use)
      self.k_sol = k_sol
      self.b_sol = b_sol


    else:
      self.set_up_bins(n_bins=n_bins)

  def generate_mock_rms_fc_list(self, n_bins = None, d_min =None,
       map_model_manager = None, model = None,
       k_sol = None,
       b_sol = None,
       n_bins_use = None,
       out = null_out()):
    mmm = map_model_manager
    if not model:
      from cctbx.development.create_models_or_maps import generate_model
      model = generate_model(n_residues = 100, b_iso = 0, log = null_out())
    # shift the model and return it with new crystal_symmetry
    from cctbx.maptbx.box import shift_and_box_model
    model = shift_and_box_model(model = model,
       crystal_symmetry = mmm.crystal_symmetry())
    model.set_shift_cart((0,0,0))
    model = mmm.remove_model_outside_map(model = model,
       boundary = 0,  return_as_new_model=True)
    map_id_model_map = 'map_for_mock_data'
    out_sav = mmm.log
    mmm.set_log(null_out())
    mmm.generate_map(model=model,
       gridding=mmm.get_any_map_manager().map_data().all(),
       d_min=d_min,
       k_sol = k_sol,
       b_sol = b_sol,
       map_id = map_id_model_map)

    self.rms_fc_list = mmm.get_rms_f_list(map_id = map_id_model_map,
       d_min = d_min,
       n_bins = n_bins).rms_f_list
    self.n_bins_use = n_bins_use  # just save it
    mmm.remove_map_manager_by_id(map_id = map_id_model_map)
    mmm.set_log(out_sav)

  def set_up_bins(self, n_bins = None):

    d_min_input_values = []
    rms_fc_values = []

    found=False
    for line in rms_fc_list.splitlines():
      if found:  # read the data
        spl=line.strip().split()
        if len(spl) != 2: continue
        d_min_input_values.append(float(spl[0]))
        rms_fc_values.append(float(spl[1]))

      elif line.strip().startswith("d_min"):
        found=True

    d_min_input_values.reverse()
    rms_fc_values.reverse()
    self.d_min_input_values=flex.double(d_min_input_values)
    self.rms_fc_values=flex.double(rms_fc_values)

    # Make a little table of values we can just use
    self.d_min = self.d_min_input_values.min_max_mean().min
    self.d_max = self.d_min_input_values.min_max_mean().max
    self.n_bins = n_bins
    self.delta_d = (self.d_max - self.d_min)/(n_bins - 1)
    self.d_min_values = flex.double()
    self.values = flex.double()
    for i in range(self.n_bins):
      self.d_min_values.append(self.d_min + i*self.delta_d)
      self.values.append(0)  # no value yet

    for d_min_1,value_1,d_min_2,value_2 in zip(
        d_min_input_values,rms_fc_values,
        d_min_input_values[1:],rms_fc_values[1:],):

      i_bin_real_1 = self.get_i_bin_real(d_min_1) # goes with value_1
      i_bin_real_2 = self.get_i_bin_real(d_min_2) # goes with value_2
      for i_bin in range(self.get_i_bin_below(d_min_1),
        self.get_i_bin_above(d_min_2)+1):
        x = (i_bin - i_bin_real_1)/max(1.e-10,i_bin_real_2 - i_bin_real_1)
        if x < -1.e-6: continue  # outside range
        if x > 1+1.e-6: continue  #
        value = value_1 + (value_2 - value_1) * x
        self.values[i_bin] = value
    assert self.values.min_max_mean().min > 0

  def maximum_d_value(self):
    return self.d_min_values[-1]

  def get_k_sol_b_sol(self):
    if not hasattr(self,'k_sol'):
      return None,None
    else:
      return self.k_sol,self.b_sol

  def get_target_scale_factors(self, f_array = None):
    if self.rms_fc_list:
      return self.rms_fc_list

    target_scale_factors = flex.double()
    for i_bin in f_array.binner().range_used():
      d_1, d_2 = f_array.binner().bin_d_range(i_bin)
      if d_1 < 0: d_1 = max(d_2 + 10, self.maximum_d_value())
      assert d_1 >= d_2
      shell_scale = self.get_expected_value_in_bin(d_2,d_1)
      target_scale_factors.append( shell_scale)
    target_scale_factors = target_scale_factors * (
        1000./max(1.e-10,target_scale_factors.min_max_mean().max))
    return target_scale_factors

  def get_expected_value_in_bin(self,d_min = None, d_max = None):
    # Return the expected value for a bin from d_min to d_max
    #  d_max = -1 for infinity
    if d_max < 0:
      d_max = self.d_max
    values = flex.double()
    for i in range(101):
      d_value = d_min + (i/100)*(d_max - d_min)
      values.append(self.get_value(d_value))
    return values.min_max_mean().mean

  def get_value(self,d_value):
    value_below = self.values[self.get_i_bin_below(d_value)]
    value_above = self.values[self.get_i_bin_above(d_value)]
    d_value_below = self.get_d_value_below(d_value)
    d_value_above = self.get_d_value_above(d_value)
    x = max(0., min(1.,(
      d_value - d_value_below)/max(1.e-10,d_value_above - d_value_below)))
    value = value_below + (value_above - value_below) * x
    return value

  def get_d_value_below(self,d_value):
    i_bin = self.get_i_bin_below(d_value)
    return self.d_min_values[i_bin]

  def get_d_value_above(self,d_value):
    i_bin = self.get_i_bin_above(d_value)
    return self.d_min_values[i_bin]

  def get_i_bin_real(self,d_value):
    return  max(0, min(self.n_bins-1,
       (self.n_bins - 1) * (d_value-self.d_min)/(self.d_max - self.d_min)))

  def get_i_bin_below(self,d_value):
    return max(0,min(self.n_bins-2, int(self.get_i_bin_real(d_value))))

  def get_i_bin_above(self,d_value):
    return max(1,min(self.n_bins-1,1 + self.get_i_bin_below(d_value)))

  def exercise(self):
    print("Exercising approx_amplitude_vs_resolution")

    for d_value in range(20):
      print ("D-min %s  rms FC  %s " %(d_value,aavr.get_value(d_value)))

    for d_min,d_max in ( (1.1,1.3), (2.7,3.0), (3.5,3.6), (10,18)):
      print ("D_min %s D_max %s  RMS FC: %s " %(
        d_min,d_max,self.get_expected_value_in_bin(d_min,d_max)))
if __name__ == "__main__":
  aavr = approx_amplitude_vs_resolution()
  aavr.exercise()

  from iotbx.map_model_manager import map_model_manager
  mmm = map_model_manager()
  mmm.generate_map()
  new_aavr = approx_amplitude_vs_resolution(d_min = 3, n_bins = 20,
    map_model_manager = mmm)
  print(list(new_aavr.get_target_scale_factors()))
