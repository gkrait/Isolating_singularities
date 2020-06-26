import flint as ft 
import sympy as sp 
import interval_arithmetic as d 
def eval_func():
  pass 
boxes=[[[1.883252480477629, 1.937869082898732], [2.768877937193937, 2.83004534774735], [0.9797068891789072, 1.003499836105708]], [[1.928471774154837, 1.96577526714105], [2.708803393800786, 2.769046902717292], [0.9664580250465434, 0.9833233824008777]], [[1.961808495228144, 1.988404527346488], [2.655966684204013, 2.708803393800787], [0.9555875038674647, 0.9679572457524716]], [[1.986374668118282, 2.007845763501218], [2.607375286616007, 2.655966684204018], [0.9461506975511593, 0.9563464330375219]], [[1.787709722387414, 1.827482810180063], [2.880646404499667, 2.904539320847356], [1.027775529940258, 1.044789677417873]], [[1.827110777969794, 1.859797576239633], [2.85737006936042, 2.887714275029508], [1.013951036250478, 1.028250198354169]], [[1.856689225381714, 1.893970741358692], [2.827912266989905, 2.859702324754173], [0.9992706008493784, 1.015273654315781]], [[2.076055217131638, 2.104378049778684], [2.3130184191659, 2.402647841591509], [0.8973226306626712, 0.9117424304002565]], [[2.10264920057907, 2.122245840467423], [2.245469856109857, 2.313018419165904], [0.8876991457645054, 0.8979392069000232]], [[2.121191402804889, 2.141431598217901], [2.17259526006482, 2.245469856109858], [0.8773299295419628, 0.8880700754126228]], [[2.140562372945778, 2.155287286825252], [2.116869675988526, 2.172595260064824], [0.869664705354579, 0.8776310657508912]], [[2.006220492213038, 2.031399220257094], [2.545214874067069, 2.607375286616012], [0.9346186363281, 0.946748672705726]], [[2.029810795421579, 2.048232102600564], [2.494775183103772, 2.545214874067074], [0.9261377823424578, 0.9352045582273846]], [[2.046517011994216, 2.078678474259009], [2.402647841591487, 2.494775183103773], [0.9107950435362318, 0.9267555105576145]], [[2.191606276974564, 2.208889379698687], [1.894297988336181, 1.964847066640629], [0.8393681490625277, 0.8491274891561583]], [[2.208303092740154, 2.221023387180586], [1.84098224604636, 1.894297988336186], [0.8322476337196034, 0.8395567854411204]], [[2.22061297516645, 2.233800276441872], [1.785196036909923, 1.840982246046361], [0.8247225122652037, 0.8323803241481474]], [[2.23344795577718, 2.244072822229732], [1.739638282564507, 1.785196036909928], [0.8185808411599355, 0.8248317814797274]], [[2.15444733765871, 2.176239628758797], [2.033257224097027, 2.116869675988527], [0.8580693288695757, 0.8699560198897954]], [[2.175233713313689, 2.192346833309379], [1.964847066640456, 2.033257224097028], [0.8488802670388866, 0.8584080472676802]], [[2.27555101546097, 2.298048876454445], [1.502584436003409, 1.597083277578766], [0.7856525486588704, 0.7992725166670633]], [[2.297097002672533, 2.31553437604673], [1.425267177787984, 1.50258443600341], [0.7745092801944568, 0.7859169230263165]], [[2.243543553543861, 2.262065045786675], [1.661233009802053, 1.739638282564511], [0.8078212337727992, 0.8187562480601492]], [[2.261398947740137, 2.276384319113106], [1.597083277578278, 1.661233009802054], [0.7990089260124775, 0.8080182079346881]], [[1.582853261362854, 1.621041415073923], [2.785191029494092, 2.821949372321393], [1.109824143856694, 1.123641506291267]], [[1.617411446023486, 1.647497019341796], [2.818765755993215, 2.850073392616441], [1.099739805823396, 1.111126809183694]], [[1.532642640420774, 1.588075855926531], [2.723875236767205, 2.785191029494093], [1.122075675931158, 1.141575478775608]], [[1.49067284682289, 1.54530686824647], [2.663109884204746, 2.732590876185757], [1.137445542543154, 1.156263485122907]], [[1.454824796164909, 1.494603520530296], [2.606261799856647, 2.663109884204747], [1.155127092966644, 1.168606076721211]], [[1.712353481996281, 1.787709722387415], [2.881389015948459, 2.918095287564693], [1.044100377845791, 1.075661439095644]], [[1.703541054498099, 1.714944531147318], [2.885252893089504, 2.8908479864046], [1.073820070289323, 1.078309918615943]], [[1.64189885010794, 1.7035410544981], [2.838158600876311, 2.896327461159873], [1.07806207364987, 1.102348174913106]], [[1.343918360969525, 1.39392703674241], [2.430260788448891, 2.506026499620245], [1.189206505514277, 1.205317364506119]], [[1.306630872342019, 1.34798005143712], [2.370339589266012, 2.432713709443261], [1.204166295113983, 1.217287890412982]], [[1.41903724928533, 1.45665817575604], [2.549621896205952, 2.606261799856648], [1.168133860513017, 1.180671065058895]], [[1.390303546220014, 1.420784659485789], [2.503265102952669, 2.55039642123811], [1.180191388040771, 1.190221663048925]], [[2.425535544901356, 2.464799519122521], [1.461819700751812, 1.520238406463993], [-0.7016549421621326, -0.6714561111846693]], [[2.393410475084039, 2.425535544901357], [1.50393645025974, 1.529903052108781], [-0.7235231698211689, -0.6999812356038038]], [[2.334747688949953, 2.409564292382313], [1.496792276019163, 1.554370342043898], [-0.766025889646805, -0.708481030803205]], [[2.45601861458602, 2.502026840636181], [1.425267177787984, 1.472099361420593], [-0.6777939169027337, -0.646033715268138]], [[2.265117493286329, 2.334747688949954], [1.464306417102085, 1.537981482804393], [-0.8086629804146261, -0.7598239890579876]], [[2.20131501964977, 2.267170019700042], [1.408381016385162, 1.494519199667703], [-0.8457424445067885, -0.8030405378563757]], [[1.835350725721682, 1.895719247694992], [1.425267177787984, 1.489497007594995], [-1.024260859610128, -0.9988507176230785]], [[1.515466066884122, 1.53820483105028], [2.22749644873323, 2.293549879534631], [-1.147719773844663, -1.139769612415749]], [[1.496757970647827, 1.515905265226848], [2.293549879534625, 2.349859882459408], [-1.154247863547448, -1.147614894372686]], [[1.553290454892071, 1.572329282300697], [2.127380911935704, 2.181667904757279], [-1.134341822666993, -1.12754979812435]], [[1.537820321207561, 1.553647901653113], [2.181667904757273, 2.227496448733231], [-1.139852667531037, -1.134250806901602]], [[1.449430929791376, 1.468506644865544], [2.435764143392512, 2.492397358623943], [-1.170476224947331, -1.164020070747089]], [[1.433398997339071, 1.449657420430106], [2.492397358623936, 2.540758240311896], [-1.175891931287928, -1.170430080697479]], [[1.468006523841825, 1.497277836015442], [2.349859882459394, 2.435764143392513], [-1.164135506097178, -1.15414639822413]], [[1.382674584494814, 1.403249476830567], [2.631775948921601, 2.692006393684931], [-1.192746556187748, -1.185995410318726]], [[1.365383363782813, 1.382930945144916], [2.692006393684925, 2.743019814985181], [-1.19840838638658, -1.192693493419236]], [[1.402804141310955, 1.433822119598953], [2.54075824031189, 2.631775948921602], [-1.186087581477825, -1.175823251118442]], [[1.308157550683802, 1.34064153659805], [2.816818468059725, 2.907016821841053], [-1.216817700321496, -1.206498997687727]], [[1.33990593451493, 1.365782589431249], [2.743019814985167, 2.816818468059726], [-1.206658834921595, -1.198326940833447]], [[1.760139349826751, 1.806848091076856], [1.555621733792544, 1.636440839762501], [-1.055653352836339, -1.036897978797245]], [[1.8006993323903, 1.846086718822156], [1.489497007594986, 1.555621733792545], [-1.039001383205987, -1.020282244540691]], [[1.696948419357231, 1.734242943207321], [1.703318019061126, 1.785056800008946], [-1.080779550989802, -1.066364772578942]], [[1.731597775497749, 1.764621355698568], [1.636440839762488, 1.703318019061127], [-1.067178147253978, -1.054128270676671]], [[1.644229762211635, 1.665801302874449], [1.86911404387335, 1.924194963766003], [-1.100993836291863, -1.092898335096372]], [[1.626856747974912, 1.644889509900654], [1.924194963765998, 1.971850674515048], [-1.10750557917578, -1.100803309381193]], [[1.664137639673702, 1.699555662211865], [1.785056800008942, 1.869114043873351], [-1.093414372375749, -1.079988870090244]], [[1.571443225024468, 1.602554155862695], [2.041839280800175, 2.127380911935705], [-1.127795588162613, -1.116600759039521]], [[1.601477466250244, 1.627545486583545], [1.971850674515035, 2.041839280800176], [-1.116873986455395, -1.107332497625868]], [[2.675995191880319, 2.709148287188679], [0.4908912468762015, 0.524349642173026], [0.4717654688810495, 0.504107439822596]], [[2.70816257444589, 2.738123499563378], [0.4754721735739187, 0.5047345442862348], [0.4447591218839698, 0.4749058658612548]], [[2.60855344525771, 2.643183346811933], [0.5316870505519107, 0.5787305840371001], [0.5303471348792494, 0.5618111368391778]], [[2.641727546078814, 2.67599519188032], [0.5093294108449989, 0.5481474711604377], [0.501704792697815, 0.533936131923425]], [[2.77724092688826, 2.809614387751462], [0.4177425180829191, 0.4570126494402861], [0.3749292563507242, 0.4063681636284753]], [[2.804668539679957, 2.834175515768741], [0.3943775761210629, 0.4213273315762874], [0.3478448415707824, 0.3780373295030455]], [[2.734838917589327, 2.778111201485537], [0.4393178303075329, 0.4981614865413235], [0.4040796125956342, 0.451758565011149]], [[2.89283659789983, 2.944365130153906], [0.1797453355001069, 0.2579713788485344], [0.2044440957303222, 0.2710569132934329]], [[2.92981790601806, 2.958404079107322], [0.1224765462756921, 0.1801665492730303], [0.1730969543942733, 0.2183997162182781]], [[2.961910299791292, 2.987197812160758], [0.04462754041896711, 0.08748545233668835], [0.1040631378174774, 0.1580928339515659]], [[2.952928766222347, 2.968501340437467], [0.08748545233668417, 0.1224765462756922], [0.1484600739916906, 0.1784060534573241]], [[2.982180264658069, 2.995772802474007], [0.009945927393446102, 0.04462754041896753], [0.0539659057668041, 0.1090743971253024]], [[2.820366672121687, 2.86921375256499], [0.3487722486757039, 0.394377576121063], [0.3090279881633347, 0.3597538176886807]], [[2.853810017021161, 2.882759508214437], [0.3102354506215276, 0.3510806342782944], [0.2896428017427739, 0.3228694787554012]], [[2.874194475687408, 2.908553977658689], [0.256433280456197, 0.3102354506215277], [0.2556309280753479, 0.2974200293917101]], [[2.970578920205049, 2.988405893904908], [0.01672222684518037, 0.04488769860680609], [-0.1403147853650912, -0.08937980116849253]], [[2.988025233109442, 2.999135417960545], [0.001450216944590311, 0.02232745445287559], [-0.08937980166533244, -0.02487423408063493]], [[2.946415986979499, 2.971774467950493], [0.03914435480842969, 0.07415999901679113], [-0.1899512911186135, -0.1403147853650895]], [[2.343990276118124, 2.365848099926291], [1.213776379377156, 1.298665937258771], [0.7417248958671206, 0.7557708869698804]], [[2.364676017210365, 2.383109921110683], [1.14432126149287, 1.213776379377157], [0.7300378953242121, 0.7421326043265716]], [[2.314793314979251, 2.331628294644757], [1.355636474829361, 1.425267177787985], [0.7641360453320051, 0.7747279114838898]], [[2.331025746889407, 2.344970708215164], [1.298665937258214, 1.355636474829362], [0.7553992355453247, 0.7643197700075292]], [[2.405901739096933, 2.422364074400928], [0.9980255675724639, 1.053780853662215], [0.7029020156583428, 0.7141136086516207]], [[2.42151910276051, 2.435555342697147], [0.9524075901767906, 0.998025567572464], [0.6935379309708474, 0.7032497017133302]], [[2.38140324404149, 2.407383998225867], [1.053780830178927, 1.144321261492873], [0.7135718206774263, 0.73074687409738]], [[2.475748113547768, 2.507027757474819], [0.7496977355251321, 0.8195204783154438], [0.6416012717777864, 0.6638583620979434]], [[2.502520400719606, 2.532653879987324], [0.6925700211603236, 0.7496977355251322], [0.621965447931285, 0.643913357759237]], [[2.433744658603749, 2.459014681095449], [0.8793196667305674, 0.9524075901767906], [0.6769513898769683, 0.6944203292900171]], [[2.456677044302328, 2.479305644402647], [0.8195204783152029, 0.8793196667305675], [0.6620291393474917, 0.678018110678026]], [[2.562060873576762, 2.589474986804852], [0.5977912433033364, 0.6305338332995262], [0.5767961893000345, 0.5978530796971427]], [[2.583637899395028, 2.612853564165159], [0.5710013857566768, 0.5977912433033366], [0.5576465706084985, 0.5804651599841776]], [[2.524738134397265, 2.569047930941044], [0.6305338293185943, 0.6925700211603252], [0.5940585890036225, 0.6264547116071282]], [[2.676596870137585, 2.697871687735659], [0.7943109409148446, 0.8694572351043965], [-0.5019107858716672, -0.4835719827482232]], [[2.660281355830595, 2.678337319894604], [0.8694572351043942, 0.9333095102109712], [-0.5163546525748776, -0.5009920022175224]], [[2.709667266138457, 2.726625935125822], [0.6805848090279195, 0.7415596734715427], [-0.4721261093953966, -0.4568676707894969]], [[2.696375099157539, 2.710807500804996], [0.7415596734715405, 0.7943109409148447], [-0.4843370976407948, -0.4715099627060307]], [[2.61948649515709, 2.637900385527291], [1.022682229502021, 1.083284081346124], [-0.5512196721530053, -0.5361706532778939]], [[2.60556012345325, 2.620745927928119], [1.08328408134612, 1.131915563663948], [-0.5628808957259103, -0.5505791196112154]], [[2.63561093677836, 2.662668519330067], [0.9333095102109592, 1.022682229502022], [-0.5373533626511828, -0.5151342878498055]], [[2.558281782280106, 2.580663484716177], [1.21763102815781, 1.276032161385307], [-0.6009923748448583, -0.5837572379784197]], [[2.541606820314024, 2.560456289604861], [1.275555749517814, 1.320689047624031], [-0.6141899675712922, -0.599821943295715]], [[2.577615581325439, 2.608338346392217], [1.131915563663944, 1.217631028157811], [-0.5852861793815742, -0.5614572790740188]], [[2.4899695723745, 2.528117112068512], [1.366504629972736, 1.425267177787985], [-0.6530803138881154, -0.6258060362362212]], [[2.521857656225531, 2.544966619695067], [1.31842618176026, 1.366504629972737], [-0.6293821282354319, -0.612144892016939]], [[2.881893539710833, 2.923002403069205], [0.1153814953755265, 0.1692598140468299], [-0.2868974558446605, -0.2357577452874414]], [[2.91621645547621, 2.949141532924421], [0.07129923458786593, 0.1156831586495081], [-0.2393362105255641, -0.1899512911186134]], [[2.838901168983457, 2.86818259872296], [0.2202386809581498, 0.2825461850204332], [-0.3402773667479331, -0.307971629948488]], [[2.862773503180642, 2.891084885455981], [0.1692598140468265, 0.2202386809581499], [-0.3120848004622018, -0.2787150712343771]], [[2.791554284331536, 2.811502004883652], [0.3758314965399165, 0.4351978859607864], [-0.3922287719909325, -0.3720951484130843]], [[2.776197910542495, 2.793259711161276], [0.4351978859607849, 0.4898468762629103], [-0.4080445102576642, -0.3912090368792646]], [[2.806282462743889, 2.845217274342616], [0.2825461850204319, 0.3758314965399166], [-0.3757695665903089, -0.3361763243351511]], [[2.737895166236414, 2.754036889750847], [0.575678945510145, 0.6322681170658899], [-0.4456726923251862, -0.4306040595256076]], [[2.725559418986216, 2.738926988604205], [0.6322681170658878, 0.6805848090279196], [-0.4574223687326879, -0.4450990075686871]], [[2.751759499898165, 2.778389735045797], [0.489846876262904, 0.5756789455101451], [-0.431950513581952, -0.4069280179947641]], [[2.150053850267501, 2.204429120572214], [1.371331180994407, 1.4364918105727], [-0.873677174252861, -0.8410122155669595]], [[2.111614269072056, 2.150356955316575], [1.353509510017523, 1.388288704293931], [-0.8937842225299646, -0.8719682679683857]], [[2.08560765334959, 2.111614269072057], [1.345658312423384, 1.361729407372393], [-0.9071444975672411, -0.893127791978631]], [[2.037325984243922, 2.086348758495839], [1.329067276289524, 1.354822052196081], [-0.9322316801831017, -0.9060936289848789]], [[1.955292488705732, 1.991659379135037], [1.341619251887513, 1.36809863179701], [-0.9713997449146541, -0.9536899190695876]], [[1.989909883383653, 2.03974296241308], [1.328895755394327, 1.353312373747569], [-0.95536847908864, -0.9299059833136238]], [[1.889332816278565, 1.955292488705733], [1.346895655597637, 1.436115434914596], [-1.002026619153407, -0.9705394684808641]], [[0.01252677538817148, 0.0603138650683445], [2.879439353229335, 2.974514737803633], [1.555717240075027, 1.567664246852382]], [[0.05985413325496745, 0.0988703373059823], [2.802122369261597, 2.87943935322934], [1.546071348527599, 1.555829768757316]], [[0.09833696178181751, 0.1382929645524683], [2.723910182702173, 2.802122369261598], [1.536199771928722, 1.546200858013405]], [[0.1377684762095851, 0.1702938081803577], [2.661074067727204, 2.723910182702178], [1.528175071624509, 1.536325735876339]], [[-0.07733652270938252, -0.0446495223120269], [3.089182370831448, 3.152770653838259], [1.581960053289938, 1.590133840064092]], [[-0.0448903929755369, -0.0130442886622301], [3.026305474680891, 3.089182370831449], [1.574057567501948, 1.582019262919605]], [[-0.01320466606618325, 0.01277988541250866], [2.974514737803624, 3.026334806561257], [1.567601389698477, 1.574097459054934]], [[1.030921864358358, 1.063565585776806], [3.499850337415542, 3.556858956728485], [-1.300609212164631, -1.291200224055728]], [[1.004659803338156, 1.032299217494783], [3.556858956724986, 3.601714823099528], [-1.30815842184547, -1.300251884650975]], [[0.9546385448312115, 1.008639005092844], [3.601714823099523, 3.675584282086779], [-1.322349361386574, -1.307120214093667]], [[1.093070723196806, 1.133424813084737], [3.356068019046911, 3.437797425449472], [-1.282503911393565, -1.270650068999962]], [[1.061970345929266, 1.094791539909514], [3.437797425449466, 3.499850337415543], [-1.291611780020637, -1.282062655372576]], [[0.8091768202666751, 0.9069884168249194], [3.73095030544818, 3.844520739044758], [-1.362941317268971, -1.335648592292274]], [[0.9003561527299625, 0.9617096717991565], [3.675584282086771, 3.74271123410899], [-1.337554256378089, -1.320493532114606]], [[0.7052971209342365, 0.768763807723133], [3.843876051671728, 3.871191815338384], [-1.39075554593726, -1.37361382300072]], [[0.6533692028473143, 0.7056689313719666], [3.860714425381281, 3.868516139609445], [-1.404500648708053, -1.390552865213586]], [[0.7681887761656623, 0.8097609096082234], [3.823469526246196, 3.85015013345232], [-1.373840153099881, -1.362513733637219]], [[0.5703808871669576, 0.6533692028473144], [3.829858137170103, 3.874752461904821], [-1.426332575307601, -1.404407376964851]], [[0.5245833406746769, 0.5706652643627674], [3.804017633583719, 3.837894567347922], [-1.43814446746416, -1.426133207407509]], [[1.25392051601231, 1.285013770532925], [2.974036399062891, 3.055949215524452], [-1.233887193672635, -1.224220220187491]], [[1.284246508697466, 1.308824545729595], [2.907016821841038, 2.974036399062892], [-1.224398168424585, -1.216663396567158]], [[1.210673267103578, 1.235441965559762], [3.106639860523367, 3.168595093323667], [-1.247246090851042, -1.239668119118755]], [[1.234898172133163, 1.254490806112254], [3.05594921552444, 3.106639860523368], [-1.239798948464645, -1.233749160689305]], [[1.166314872062655, 1.191798285546401], [3.217346659907928, 3.27693190792774], [-1.260719750849207, -1.25304615274301]], [[1.191165321640507, 1.211176339055835], [3.168595093323661, 3.217346659907929], [-1.253202341545049, -1.247122354992473]], [[1.131487219138907, 1.167603384764931], [3.276931907892302, 3.356068019046912], [-1.271141532120167, -1.260397692615998]], [[0.450731917395545, 0.524583340674677], [3.737236857204096, 3.813807940647724], [-1.457213170565211, -1.438089288420697]], [[0.3997106387733055, 0.4522392959773724], [3.687122356329893, 3.740388954474258], [-1.470185650701914, -1.456786254894892]], [[0.3432975394884848, 0.4065946092635829], [3.614798890028737, 3.689478227901574], [-1.48453493554697, -1.46847443890771]], [[0.3034457106805704, 0.3476957851711536], [3.553869075201785, 3.614798890028738], [-1.494636801241454, -1.483445587815372]], [[0.2461895160828349, 0.3068321475967375], [3.464447248291217, 3.55522072617614], [-1.50908521619486, -1.493808951474127]], [[0.2038752467231475, 0.248915314346961], [3.39015168013102, 3.464447248291218], [-1.519735544934522, -1.508416759032594]], [[0.1584673959106966, 0.2057180029526086], [3.308445829337261, 3.390931955794888], [-1.531135010223486, -1.519285182988713]], [[0.1226052588287029, 0.1595583614564186], [3.240705576477624, 3.308445829337262], [-1.540124408660009, -1.530868747495722]], [[0.08232521064418975, 0.1232207559953372], [3.163441764381447, 3.240705576477625], [-1.550208040716588, -1.539975487771074]], [[0.05028308447339331, 0.0828303064668111], [3.100313057402689, 3.163441764381448], [-1.558223769517915, -1.550084391889261]], [[0.01763254564310245, 0.0506046302040621], [3.035488330844447, 3.100440255244453], [-1.566387936810696, -1.558144595311572]], [[-0.008898472105635389, 0.01780001837088771], [2.982307090699998, 3.035488330844448], [-1.573020910866189, -1.566346354320475]], [[-0.04159139045910103, -0.00876585684234843], [2.91678456191558, 2.982307090700003], [-1.581194380104902, -1.572987917874853]], [[-0.06827021751182795, -0.04141256897441997], [2.86317522018287, 2.916784561915581], [-1.587866319492534, -1.581150467763988]], [[0.4618648395935464, 0.5393954421461323], [2.077004365732816, 2.158787653234855], [1.434349780653566, 1.45426201880633]], [[0.5305916487471077, 0.5926027551871091], [2.022986243443111, 2.07700436573282], [1.420472285442062, 1.436539250789069]], [[0.6493537422639726, 0.7233325376530181], [1.927009208030491, 1.979627997747961], [1.385815510869417, 1.405614398456958]], [[0.588827879290803, 0.6493537422639727], [1.969731820288429, 2.026411794025806], [1.405487160457522, 1.42146533607385]], [[0.7713553911447023, 0.8305870623164261], [1.913118751396499, 1.925286052431467], [1.356811914724402, 1.373028262767621]], [[0.7229231727166625, 0.7713553911447024], [1.918713610734013, 1.938198240025033], [1.372926873516343, 1.386013763701811]], [[0.2064733929963208, 0.2535523021832424], [2.503927149190693, 2.589352550746268], [1.507249464081779, 1.51908220204836]], [[0.1696324467460728, 0.2076552072002623], [2.589352550746259, 2.661074067727205], [1.518795253952585, 1.528335242298389]], [[0.2826474273760479, 0.3205898243203366], [2.384816539756408, 2.448970426742682], [1.490313210142191, 1.499894423324984]], [[0.2525925030428078, 0.2836786586767603], [2.448970426742674, 2.503927149190697], [1.499643501153973, 1.507477488850724]], [[0.355146692589751, 0.4026517772199786], [2.253093541312346, 2.323794935283304], [1.469465368798162, 1.48153166954188]], [[0.3192800838202242, 0.3572353596092963], [2.323794935283296, 2.384816539756409], [1.481019540129733, 1.490630738266158]], [[0.3971615256352184, 0.470804620616049], [2.158787653234849, 2.25309354131235], [1.452046389548871, 1.470826671340181]], [[1.209405981462292, 1.254131799153143], [2.224400938900139, 2.285127577443699], [1.233938581153283, 1.24760747644072]], [[1.178713517812621, 1.213835986363792], [2.181414611025879, 2.227066260416622], [1.246345232817914, 1.25696462124741]], [[1.246810833399831, 1.312174056522173], [2.280804268358795, 2.372229769883475], [1.215786935362243, 1.236028811915719]], [[1.033851450340963, 1.13610854495627], [2.026169121780117, 2.10818222458316], [1.270171766298136, 1.299636495422085]], [[1.113951884519775, 1.185511132643453], [2.10375986878587, 2.18141461102588], [1.255140076589178, 1.276287576211084]], [[0.9634124616126288, 1.044269958469339], [1.95303812551386, 2.032197553996884], [1.296677955683544, 1.320126707400716]], [[0.8299733174528409, 0.9047551835271502], [1.909242591387182, 1.939810929414247], [1.336311467185174, 1.357160519999176]], [[0.9040184122241364, 0.9642895676051298], [1.928841518512021, 1.968028608150939], [1.319613809450584, 1.336682638094038]], [[-0.1170006506397822, -0.07675397619811788], [3.152553085115033, 3.228937326054354], [1.589990592150522, 1.600059673873985]], [[-0.1170669211773508, -0.06776245890122387], [2.766129726669345, 2.863175220182881], [-1.600074822626443, -1.587741629880258]], [[2.997818933471124, 2.999684738993313], [0.0006592781692382269, 0.004797233458152947], [0.014546, 0.03813463518828137]], [[2.995634816897783, 2.997755709752326], [0.004961849223200528, 0.009945964086264336], [0.03875591217760954, 0.05396600000000001]], [[2.997744308662616, 2.997821151022688], [0.004797233458151311, 0.004973473671868386], [0.03813463518827553, 0.03879936256298878]], [[2.999682626422045, 3.000000044762371], [-0.0, 0.0006592781692382679], [-0.005164000000000001, 0.01454600000000001]], [[2.999535905014029, 2.999960095455392], [7.625983501450428e-05, 0.0009153424085324696], [-0.01758991497585444, -0.005164]], [[2.999072088132142, 2.999493589984398], [0.0009467162785655198, 0.001765128797156666], [-0.02487400000000001, -0.01838127438250336]], [[2.999493303093181, 2.999536018416645], [0.0008833406806067206, 0.000964582463792688], [-0.01838127438250404, -0.0175899149758418]]]
ftboxes=[ [d.ftconstructor(Bi[0],Bi[1]) for Bi in B ] for B in boxes ] 
n=len(boxes[0])
m=len(boxes)
m1=[]
m2=[]
for B in ftboxes: 
   m1.append((-ft.arb.sin(B[2])**4 - 3)*ft.arb.sin(B[2]) + 4*ft.arb.sin(B[2])**3*ft.arb.cos(B[2])**2) 
   m2.append(-2*(ft.arb.sin(8*B[2]) + 3)*ft.arb.sin(B[2])*ft.arb.cos(B[2]) - 8*ft.arb.sin(B[2])**2*ft.arb.cos(8*B[2])) 
innrer_loops=[i for i in range(m) if 0 in m1[i] and 0 in m2[i] ]
print(innrer_loops)
