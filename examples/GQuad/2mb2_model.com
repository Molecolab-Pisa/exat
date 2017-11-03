%nosave
#p td=nstates=10 m062x 6-31G(d) eet=(fragment=12,fragmentcavity) scrf=(iefpcm,solvent=water,read) nosymm integral=ultrafinegrid

2mb2 model1 G projected

0 1
C(fragment=1)      2.876099     8.129722     0.726309
H(fragment=1)      3.312598     7.900521    -0.257746
H(fragment=1)      3.685127     8.415887     1.415520
H(fragment=1)      2.165461     8.962816     0.627930
N(fragment=1)      2.164522     6.972809     1.241060
C(fragment=1)      0.827240     6.869354     1.570651
H(fragment=1)      0.153438     7.720319     1.458567
N(fragment=1)      0.478255     5.665787     2.019046
C(fragment=1)      1.660233     4.953635     1.975238
C(fragment=1)      1.932347     3.578009     2.338349
O(fragment=1)      1.207063     2.695524     2.771628
N(fragment=1)      3.325685     3.306351     2.111812
H(fragment=1)      3.587621     2.350516     2.347893
C(fragment=1)      4.279612     4.175082     1.631865
N(fragment=1)      5.548859     3.683970     1.500547
H(fragment=1)      5.792075     2.734677     1.738781
H(fragment=1)      6.257215     4.313040     1.149089
N(fragment=1)      4.037007     5.427687     1.300075
C(fragment=1)      2.719177     5.741792     1.496727
C(fragment=2)      5.210667     4.968654    -2.716202
H(fragment=2)      5.520492     4.459492    -3.641508
H(fragment=2)      6.022651     4.887726    -1.977581
H(fragment=2)      5.012050     6.028820    -2.929109
N(fragment=2)      3.996366     4.365995    -2.194634
C(fragment=2)      2.770345     4.964203    -1.978588
H(fragment=2)      2.613160     6.020447    -2.202989
N(fragment=2)      1.850370     4.136905    -1.488031
C(fragment=2)      2.524193     2.936503    -1.382446
C(fragment=2)      2.061392     1.647393    -0.911007
O(fragment=2)      0.975017     1.281127    -0.488885
N(fragment=2)      3.145294     0.706953    -0.996085
H(fragment=2)      2.887647    -0.226070    -0.678203
C(fragment=2)      4.425253     0.944660    -1.443543
N(fragment=2)      5.288367    -0.115671    -1.431133
H(fragment=2)      5.017633    -1.033643    -1.113176
H(fragment=2)      6.230653     0.047450    -1.757489
N(fragment=2)      4.851206     2.115837    -1.873510
C(fragment=2)      3.855198     3.053095    -1.814139
C(fragment=3)      6.470246     0.446363    -5.669412
H(fragment=3)      6.650368    -0.317575    -6.441023
H(fragment=3)      7.040610     0.180234    -4.766479
H(fragment=3)      6.803242     1.426972    -6.038530
N(fragment=3)      5.052736     0.524483    -5.362293
C(fragment=3)      4.209767     1.611233    -5.488914
H(fragment=3)      4.579193     2.568726    -5.859650
N(fragment=3)      2.955126     1.360781    -5.122073
C(fragment=3)      2.993532     0.035480    -4.737038
C(fragment=3)      1.936847    -0.819623    -4.236579
O(fragment=3)      0.752957    -0.605317    -4.024389
N(fragment=3)      2.481469    -2.124165    -3.975685
H(fragment=3)      1.788218    -2.780090    -3.619133
C(fragment=3)      3.781837    -2.538551    -4.155539
N(fragment=3)      4.063696    -3.837504    -3.835071
H(fragment=3)      3.365671    -4.474590    -3.483008
H(fragment=3)      5.015218    -4.152487    -3.963180
N(fragment=3)      4.746934    -1.764306    -4.610617
C(fragment=3)      4.282099    -0.504511    -4.876802
C(fragment=4)      8.640573    -1.853049     1.468410
H(fragment=4)      8.395356    -2.626403     0.724622
H(fragment=4)      8.927889    -2.344810     2.410242
H(fragment=4)      9.478767    -1.243359     1.101749
N(fragment=4)      7.496544    -0.984633     1.684765
C(fragment=4)      7.409214     0.379161     1.484526
H(fragment=4)      8.264776     0.949338     1.118979
N(fragment=4)      6.214450     0.887474     1.777068
C(fragment=4)      5.491185    -0.213843     2.189548
C(fragment=4)      4.117732    -0.310431     2.639594
O(fragment=4)      3.247220     0.536278     2.773165
N(fragment=4)      3.830755    -1.682240     2.959427
H(fragment=4)      2.875574    -1.822928     3.284697
C(fragment=4)      4.684878    -2.758051     2.869306
N(fragment=4)      4.180764    -3.976126     3.232056
H(fragment=4)      3.232319    -4.098757     3.552182
H(fragment=4)      4.799049    -4.772913     3.169831
N(fragment=4)      5.935444    -2.675326     2.460399
C(fragment=4)      6.263814    -1.385264     2.141146
C(fragment=5)      5.270493    -6.244191    -0.377189
H(fragment=5)      4.723470    -6.949960    -1.020682
H(fragment=5)      5.302263    -6.648275     0.646107
H(fragment=5)      6.295409    -6.120829    -0.755420
N(fragment=5)      4.615403    -4.947824    -0.384348
C(fragment=5)      5.132897    -3.731295    -0.784237
H(fragment=5)      6.156244    -3.650292    -1.154257
N(fragment=5)      4.274410    -2.720775    -0.669029
C(fragment=5)      3.138433    -3.322676    -0.165618
C(fragment=5)      1.853729    -2.747520     0.176394
O(fragment=5)      1.436575    -1.600958     0.113726
N(fragment=5)      0.994121    -3.791533     0.663858
H(fragment=5)      0.068507    -3.454855     0.924340
C(fragment=5)      1.296245    -5.127570     0.800765
N(fragment=5)      0.307450    -5.935085     1.290283
H(fragment=5)     -0.604173    -5.586084     1.543676
H(fragment=5)      0.518267    -6.917747     1.394936
N(fragment=5)      2.463719    -5.655762     0.490961
C(fragment=5)      3.325017    -4.702002     0.019246
C(fragment=6)      0.644876    -8.123838    -2.145142
H(fragment=6)     -0.142795    -8.618534    -2.733597
H(fragment=6)      0.425323    -8.257575    -1.074937
H(fragment=6)      1.618461    -8.577436    -2.379899
N(fragment=6)      0.702142    -6.709859    -2.472468
C(fragment=6)      1.764075    -5.995446    -2.991668
H(fragment=6)      2.715331    -6.482817    -3.211427
N(fragment=6)      1.498980    -4.705080    -3.181678
C(fragment=6)      0.190037    -4.582547    -2.760009
C(fragment=6)     -0.669710    -3.417470    -2.719076
O(fragment=6)     -0.474599    -2.254648    -3.038959
N(fragment=6)     -1.951027    -3.805807    -2.195694
H(fragment=6)     -2.608866    -3.029906    -2.139226
C(fragment=6)     -2.342192    -5.061275    -1.788499
N(fragment=6)     -3.621680    -5.186422    -1.323018
H(fragment=6)     -4.260972    -4.408016    -1.273448
H(fragment=6)     -3.919531    -6.103660    -1.021390
N(fragment=6)     -1.563625    -6.124549    -1.823780
C(fragment=6)     -0.325400    -5.810896    -2.316150
C(fragment=7)     -0.920693    -7.047968     5.313080
H(fragment=7)     -1.776480    -7.242600     4.648870
H(fragment=7)     -1.295724    -6.871171     6.332650
H(fragment=7)     -0.250760    -7.919701     5.314097
N(fragment=7)     -0.175544    -5.892216     4.845373
C(fragment=7)      1.147525    -5.829973     4.453878
H(fragment=7)      1.780453    -6.718726     4.469215
N(fragment=7)      1.536725    -4.614985     4.074736
C(fragment=7)      0.397308    -3.851259     4.231592
C(fragment=7)      0.180226    -2.440922     3.982464
O(fragment=7)      0.928652    -1.565212     3.575598
N(fragment=7)     -1.184139    -2.119170     4.301045
H(fragment=7)     -1.406852    -1.136865     4.147945
C(fragment=7)     -2.157604    -2.973486     4.767446
N(fragment=7)     -3.391362    -2.432151     5.000296
H(fragment=7)     -3.595852    -1.457115     4.843597
H(fragment=7)     -4.113628    -3.050215     5.342859
N(fragment=7)     -1.965039    -4.257533     4.995715
C(fragment=7)     -0.676747    -4.620093     4.707857
C(fragment=8)     -5.434390    -5.099909     2.736344
H(fragment=8)     -6.237501    -4.939996     2.000907
H(fragment=8)     -5.728554    -4.635352     3.689802
H(fragment=8)     -5.281741    -6.178500     2.884900
N(fragment=8)     -4.192574    -4.518437     2.257238
C(fragment=8)     -2.991330    -5.153616     2.009818
H(fragment=8)     -2.879985    -6.227246     2.169926
N(fragment=8)     -2.034363    -4.337583     1.574240
C(fragment=8)     -2.656605    -3.105609     1.539528
C(fragment=8)     -2.137175    -1.810897     1.149457
O(fragment=8)     -1.033715    -1.465818     0.754684
N(fragment=8)     -3.181246    -0.832458     1.287595
H(fragment=8)     -2.882693     0.106465     1.028263
C(fragment=8)     -4.472892    -1.042700     1.714755
N(fragment=8)     -5.290725     0.051911     1.763762
H(fragment=8)     -4.979725     0.975276     1.503485
H(fragment=8)     -6.241097    -0.090838     2.075905
N(fragment=8)     -4.950300    -2.218836     2.070774
C(fragment=8)     -3.994076    -3.191955     1.958149
C(fragment=9)     -8.523758    -1.743326    -0.321083
H(fragment=9)     -9.208936    -1.279309    -1.046787
H(fragment=9)     -8.458702    -1.098217     0.568316
H(fragment=9)     -8.909766    -2.731049    -0.031029
N(fragment=9)     -7.208479    -1.915969    -0.912687
C(fragment=9)     -6.506522    -3.089461    -1.107178
H(fragment=9)     -6.924705    -4.051019    -0.805113
N(fragment=9)     -5.317722    -2.913292    -1.679142
C(fragment=9)     -5.251329    -1.546596    -1.863222
C(fragment=9)     -4.200427    -0.735150    -2.442362
O(fragment=9)     -3.113412    -1.032776    -2.913921
N(fragment=9)     -4.607493     0.643114    -2.407811
H(fragment=9)     -3.910716     1.273183    -2.802118
C(fragment=9)     -5.787650     1.154269    -1.916927
N(fragment=9)     -5.950219     2.509684    -1.993387
H(fragment=9)     -5.250444     3.120051    -2.386878
H(fragment=9)     -6.812872     2.895821    -1.635970
N(fragment=9)     -6.746946     0.420095    -1.389004
C(fragment=9)     -6.410801    -0.906917    -1.396357
C(fragment=10)     -6.643886     2.905658     4.729478
H(fragment=10)     -6.920103     3.367231     3.769227
H(fragment=10)     -6.486148     3.701337     5.473415
H(fragment=10)     -7.454891     2.243891     5.065732
N(fragment=10)     -5.434929     2.115208     4.576598
C(fragment=10)     -5.264790     0.762622     4.798448
H(fragment=10)     -6.095729     0.139077     5.132285
N(fragment=10)     -4.026138     0.334118     4.567046
C(fragment=10)     -3.360848     1.478041     4.173370
C(fragment=10)     -1.977799     1.663617     3.784808
O(fragment=10)     -1.046727     0.876987     3.701939
N(fragment=10)     -1.769917     3.048545     3.460358
H(fragment=10)     -0.812473     3.250823     3.176898
C(fragment=10)     -2.698163     4.064416     3.497565
N(fragment=10)     -2.261486     5.310362     3.142101
H(fragment=10)     -1.309770     5.494219     2.863711
H(fragment=10)     -2.934509     6.063687     3.165790
N(fragment=10)     -3.957780     3.900644     3.850421
C(fragment=10)     -4.212423     2.594438     4.171346
C(fragment=11)     -5.444369     6.019718     0.159807
H(fragment=11)     -5.393350     6.367102    -0.883337
H(fragment=11)     -4.985278     6.779966     0.809950
H(fragment=11)     -6.495532     5.876054     0.448159
N(fragment=11)     -4.748435     4.752638     0.301035
C(fragment=11)     -5.258980     3.539547     0.719853
H(fragment=11)     -6.309271     3.437772     0.997782
N(fragment=11)     -4.359015     2.559244     0.741988
C(fragment=11)     -3.202088     3.178130     0.312277
C(fragment=11)     -1.871712     2.637078     0.123585
O(fragment=11)     -1.422931     1.512974     0.289515
N(fragment=11)     -1.006552     3.686964    -0.340900
H(fragment=11)     -0.049392     3.373936    -0.495246
C(fragment=11)     -1.340432     5.000740    -0.580767
N(fragment=11)     -0.338365     5.819676    -1.021935
H(fragment=11)      0.604632     5.494174    -1.169826
H(fragment=11)     -0.572223     6.786206    -1.201239
N(fragment=11)     -2.549507     5.497874    -0.410600
C(fragment=11)     -3.417984     4.537134     0.032949
C(fragment=12)     -2.661453     6.708834    -4.217280
H(fragment=12)     -2.330707     6.856331    -5.256645
H(fragment=12)     -1.940887     7.196161    -3.542948
H(fragment=12)     -3.655070     7.159124    -4.080552
N(fragment=12)     -2.750564     5.290761    -3.915803
C(fragment=12)     -3.860084     4.569468    -3.520383
H(fragment=12)     -4.829972     5.053727    -3.395789
N(fragment=12)     -3.614372     3.276845    -3.319947
C(fragment=12)     -2.268125     3.160077    -3.603027
C(fragment=12)     -1.406962     1.995853    -3.568051
O(fragment=12)     -1.632774     0.828975    -3.285290
N(fragment=12)     -0.078309     2.391209    -3.948684
H(fragment=12)      0.583074     1.616272    -3.945459
C(fragment=12)      0.351507     3.651897    -4.296452
N(fragment=12)      1.672674     3.783343    -4.623048
H(fragment=12)      2.314890     3.005808    -4.614949
H(fragment=12)      1.999131     4.704450    -4.879966
N(fragment=12)     -0.428167     4.714422    -4.329798
C(fragment=12)     -1.710825     4.394152    -3.974517

alpha=1.8
