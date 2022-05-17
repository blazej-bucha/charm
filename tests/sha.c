/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






/* Function prototypes */
/* ------------------------------------------------------------------------- */
int cmp_arrays(REAL *, REAL *, size_t, REAL);
/* ------------------------------------------------------------------------- */






#define ncoeffs 66






/* Function to test the "sha" module */
/* ------------------------------------------------------------------------- */
int sha(unsigned long nmax, char SHCs_file[])
{
    /* --------------------------------------------------------------------- */
    if (nmax != 10)
    {
        fprintf(stderr, "\"nmax\" has to be \"10\".\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Initialize an error structure to be used throughout the function */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize an \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Initialize a "shc" structure and read test coefficients from a text file
     * */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_ref = CHARM(shc_init)(nmax, ADDP(1.0), ADDP(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }
    FILE *fptr = fopen(SHCs_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n", SHCs_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_mtx)(fptr, nmax, shcs_ref, err);
    CHARM(err_handler)(err, 1);
    fclose(fptr);
    /* --------------------------------------------------------------------- */






    /* Check spherical harmonic analysis with points values using all supported
     * quadratures */
    /* --------------------------------------------------------------------- */
    CHARM(crd) *grd;
    REAL *f;
    CHARM(shc) *shcs_out;
    REAL *dda = (REAL *)malloc((nmax + 1) * sizeof(REAL));
    if (dda == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances.\n");
        exit(CHARM_FAILURE);
    }
    REAL *dda_ref = (REAL *)calloc((nmax + 1), sizeof(REAL));
    if (dda_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesize in "delta_r" height above the reference sphere, then
     * do the analysis on the same sphere and rescale the coefficients to
     * "shcs_ref->r", so that the results can be compared with respect to the
     * input coefficients. */
    REAL delta_r = ADDP(1000.0);


    /* Loop over all supported quadratures */
    int errnum = 0;
    for (int i = 0; i < 3; i++)
    {
        if (i == 0)
            printf("    Surface SHA with a Gauss--Legendre point grid...\n");
        else if (i == 1)
            printf("    Surface SHA with a Driscoll--Healy point grid "
                   "(DH1)...\n");
        else if (i == 2)
            printf("    Surface SHA with a Driscoll--Healy point grid "
                   "(DH2)...\n");


        /* Now loop over maximum harmonic degrees */
        for (unsigned long nmax_tmp = 0; nmax_tmp <= nmax; nmax_tmp++)
        {
            /* Create the grids on a sphere with the radius "shcs_ref->r
             * + delta_r" */
            if (i == 0)
            {
                grd = CHARM(crd_gl)(nmax_tmp, shcs_ref->r + delta_r);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Gauss--Legendre "
                                    "grid.\n");
                    exit(CHARM_FAILURE);
                }
            }
            else if (i == 1)
            {
                grd = CHARM(crd_dh1)(nmax_tmp, shcs_ref->r + delta_r);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Driscoll--Healy "
                                    "grid (DH1).\n");
                    exit(CHARM_FAILURE);
                }
            }
            else if (i == 2)
            {
                grd = CHARM(crd_dh2)(nmax_tmp, shcs_ref->r + delta_r);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Driscoll--Healy "
                                    "grid (DH2).\n");
                    exit(CHARM_FAILURE);
                }
            }


            shcs_out = CHARM(shc_init)(nmax_tmp, shcs_ref->mu, shcs_ref->r);
            if (shcs_out == NULL)
            {
                fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
                exit(CHARM_FAILURE);
            }


            f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
            if (f == NULL)
            {
                fprintf(stderr, "Failed to initialize an array to store the "
                                "synthesized signal.\n");
                exit(CHARM_FAILURE);
            }


            /* We do the synthesis on the sphere with the radius "shcs_ref->r
             * + delta_r" */
            CHARM(shs_point)(grd, shcs_ref, nmax_tmp, f, err);
            CHARM(err_handler)(err, 1);


            /* Now do the analysis on the sphere with the radius "shcs_out->r
             * + delta_r" and rescale the coefficients to "shcs_out->r ==
             * shcs_ref->r", so that the coefficients can be validated with
             * respect to "shcs_ref" */
            CHARM(sha_point)(grd, f, nmax_tmp, shcs_out, err);
            CHARM(err_handler)(err, 1);


            CHARM(shc_dda)(shcs_ref, shcs_out, nmax_tmp, dda, err);
            CHARM(err_handler)(err, 1);


            errnum += cmp_arrays(dda, dda_ref, nmax_tmp + 1,
                                 CHARM(glob_threshold2));


            CHARM(shc_free)(shcs_out);
            CHARM(crd_free)(grd);
            free(f);
        }
    }
    /* --------------------------------------------------------------------- */






    /* Check spherical harmonic analysis with block-mean values in cells using
     * all supported quadratures */
    /* --------------------------------------------------------------------- */
    /* In these tests, we remove a few low-degree coefficients, so that we do
     * not loose numerical precision (the coefficient is several order of
     * magnitude larger than the other coefficients).  It is OK to have to
     * coefficients in quadruple precision, but it is not ideal in single and
     * double precision. */
    /* ..................................................................... */
    shcs_ref->c[0][0] = ADDP(0.0);
    shcs_ref->c[2][0] = ADDP(0.0);
    /* ..................................................................... */


    /* Symmetric grid */
    /* ..................................................................... */
    grd = CHARM(crd_init)(CHARM_CRD_CELLS_GRID, 15, 30);
    if (grd == NULL)
    {
        fprintf(stderr, "Failed to initialize the grid cells.\n");
        exit(CHARM_FAILURE);
    }


    f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
    if (f == NULL)
    {
        fprintf(stderr, "Failed to initialize an array to store the "
                        "synthesized signal.\n");
        exit(CHARM_FAILURE);
    }




    /* Define grid cells */
    for (size_t i = 0; i < grd->nlat; i++)
    {
        grd->lat[2 * i]      = -PI_2 + (REAL)i       * PI / grd->nlat;
        grd->lat[2 * i + 1]  = -PI_2 + (REAL)(i + 1) * PI / grd->nlat;
        grd->r[i] = shcs_ref->r + delta_r;
    }
    for (size_t j = 0; j < grd->nlon; j++)
    {
        grd->lon[2 * j]      = (REAL)j       * (ADDP(2.0) * PI) / grd->nlon;
        grd->lon[2 * j + 1]  = (REAL)(j + 1) * (ADDP(2.0) * PI) / grd->nlon;
    }


    /* Do the synthesis in the cells to get a signal to be harmonically
     * analyzed */
    CHARM(shs_cell)(grd, shcs_ref, nmax, f, err);
    CHARM(err_handler)(err, 1);


    shcs_out = CHARM(shc_init)(nmax, shcs_ref->mu, shcs_ref->r);
    if (shcs_out == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(sha_cell)(grd, f, nmax, CHARM_SHA_CELL_AQ, shcs_out, err);
    CHARM(err_handler)(err, 1);


    REAL shcs_out_symm_ref_cnm[ncoeffs] = {
                      shcs_out->c[0][0], /* We do not want to validate this
                                            coefficient, since it has been
                                            removed from the input
                                            coefficients.  It is non-zero only
                                            due to numerical errors and is
                                            a bit difficult to validate
                                            properly.  By taking the value from
                                            the recovered set of coefficients,
                                            we make sure that it always passes
                                            the comparison test. */
                      ADDP( 4.7454677094139766719734358859253947e-11),
                      ADDP(-4.7369019130196063124846136229108182e-04),
                      ADDP( 9.1639562352415356457832071949805106e-07),
                      ADDP( 4.1995801383664354376264320003975419e-07),
                      ADDP( 6.2085971552317452190352330050077898e-08),
                      ADDP(-2.3249697099775248635905878380426906e-07),
                      ADDP( 7.4951727837435956514459913331502198e-08),
                      ADDP(-9.1774493939302794859530196484006215e-08),
                      ADDP( 2.1198558488022290034957280477261947e-08),
                      ADDP(-1.2682812077535012035626930299710638e-07),
                      ADDP( 7.9475652657652941955556203439722187e-09),
                      ADDP(-3.8637271217486361448925654890238139e-09),
                      ADDP( 1.9648188607985963699971331635341456e-06),
                      ADDP(-5.0530703279271140043469359199720947e-07),
                      ADDP(-2.8118575261760469260641559512334829e-08),
                      ADDP(-7.2275198843562173535315577390005510e-08),
                      ADDP( 2.5457980883798820101580461464929723e-07),
                      ADDP( 1.2400315486844960968477194707361625e-08),
                      ADDP( 1.2803284688615461021827459179449718e-07),
                      ADDP( 5.0867896357477706041090293469668613e-08),
                      ADDP( 2.3655499657831677664284854663688666e-09),
                      ADDP( 8.8563870045644499469911299577819729e-07),
                      ADDP( 3.3306533720793623753439448136709878e-07),
                      ADDP( 6.1886447807467333759124036795222872e-07),
                      ADDP( 4.9220785897406328989076514911835370e-08),
                      ADDP( 3.0259724846128266618554084233472964e-07),
                      ADDP( 6.6470351595466854576299293624550974e-08),
                      ADDP( 4.3664012880081215833686829719560939e-08),
                      ADDP(-5.9918381460448039126844980078535215e-08),
                      ADDP( 6.9057754465292816559257742546319902e-07),
                      ADDP( 9.3698987516828215714089972528303072e-07),
                      ADDP(-4.0624305004999227658905543986658771e-07),
                      ADDP( 6.6167554216667947265459593486801383e-08),
                      ADDP( 2.0025388296889353564801957820184471e-07),
                      ADDP(-6.4150249914332240435572178624414256e-09),
                      ADDP(-1.1894775229298843533632281855352478e-07),
                      ADDP(-3.4233756503868877836807717465844912e-10),
                      ADDP(-1.7851776402172440097647171807826152e-07),
                      ADDP(-2.7599921227433063245914241363228722e-07),
                      ADDP(-8.6960739690732369915965065668061294e-08),
                      ADDP(-2.4224119038455644012993687608648621e-07),
                      ADDP(-2.0803261266141500877467050696504429e-07),
                      ADDP(-1.9890770297760397003258982476623938e-08),
                      ADDP(-7.4791393812826393182248042482862931e-08),
                      ADDP( 1.5808199854881863743594694820264574e-07),
                      ADDP(-2.3634105414120480594824104194818242e-07),
                      ADDP( 3.3422400163097533262153683787679425e-09),
                      ADDP(-2.9261877985108827339843653371081135e-08),
                      ADDP(-1.2256436067587534470802015972773190e-08),
                      ADDP(-3.9809472814960087218447636025810532e-08),
                      ADDP( 7.2828760291976033637040760490193230e-09),
                      ADDP(-3.0015663157343560954801109719422561e-07),
                      ADDP(-5.4571392397950643020449016974203620e-08),
                      ADDP( 3.9394218015314553586501074369600597e-08),
                      ADDP(-3.0154898141819735574021476299494347e-08),
                      ADDP(-9.7794735528802857644586797865389838e-10),
                      ADDP( 5.3808134060485507990707698201113362e-08),
                      ADDP(-9.0568876444962032778684416263081202e-08),
                      ADDP( 7.8136383624837493740876984596636453e-09),
                      ADDP(-9.5565684356907844255008338788979383e-08),
                      ADDP( 1.4095428243382017333713512686424309e-07),
                      ADDP( 2.6908437351803154057502941325105060e-08),
                      ADDP(-3.4617231922453544280917180314029827e-08),
                      ADDP( 8.7429016507122066708469112566959367e-08),
                      ADDP( 6.7535689277504987473850868046495163e-08)};
    REAL shcs_out_symm_ref_snm[ncoeffs] = {
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 0.0000000000000000000000000000000000e+00),
                      ADDP( 8.8492304573374505049632524211856023e-10),
                      ADDP(-2.2473189538218186423838318793422686e-09),
                      ADDP( 2.4020779004998599792470171108364742e-07),
                      ADDP(-4.4646176109047101681109795665388961e-07),
                      ADDP(-8.1500166244169440938761828435820084e-08),
                      ADDP( 1.5210197225909814421348020270553901e-08),
                      ADDP( 8.0846459970677777428733352649620257e-08),
                      ADDP( 3.7969407894150146672099047513103616e-08),
                      ADDP( 1.8845554574026741954927429264833130e-08),
                      ADDP(-9.3375024287163768645160426840026865e-08),
                      ADDP(-1.3719578944287833362700654125094620e-06),
                      ADDP(-6.0236088591312079486306028079522953e-07),
                      ADDP( 6.1259741202868398544369888911055984e-07),
                      ADDP(-3.0306973953100073203110032992872384e-07),
                      ADDP(-3.2391374716595540991809053470112229e-07),
                      ADDP( 6.5850887865819533754156123606696507e-08),
                      ADDP( 4.4304629066850852529047476631006797e-08),
                      ADDP(-2.9710979000573484226163643421544768e-08),
                      ADDP(-3.8797060855931837941392464291092109e-08),
                      ADDP( 1.3582158417279769588526621453658476e-06),
                      ADDP(-1.9110242505282179358153346898388849e-07),
                      ADDP(-1.9175804928134839864837128004250500e-07),
                      ADDP(-6.5168043906141878095949239431915418e-10),
                      ADDP(-1.8875258360928175050358212810752111e-07),
                      ADDP(-7.8815021583245248667596106444556623e-08),
                      ADDP(-6.5669935549555061015268700025751724e-08),
                      ADDP(-1.1626005751814777410890843570432805e-07),
                      ADDP( 2.8410819765910957010790280706434485e-07),
                      ADDP( 4.3188358604635117692891636139397085e-08),
                      ADDP(-4.1506317787013936406890489212754024e-07),
                      ADDP(-1.0448652771806180068491531385389123e-07),
                      ADDP( 4.1302801037011390484683196221424121e-08),
                      ADDP( 1.0882960061065380342335103646975707e-08),
                      ADDP(-6.0794272356712211749764422849153770e-08),
                      ADDP(-6.0489330664849967973649692769822587e-07),
                      ADDP(-4.7100298483654507557201412494288839e-07),
                      ADDP( 4.7188269982450173692302627221577240e-09),
                      ADDP( 5.8417295813122735535644449123810345e-08),
                      ADDP(-4.4244438926950942596880666217527800e-08),
                      ADDP(-3.7883074157345382260136647370051960e-08),
                      ADDP(-2.0115579073486506478436007416285296e-07),
                      ADDP( 1.3341433494942722794652206704345337e-07),
                      ADDP( 2.4440340067449308539502769126431027e-07),
                      ADDP( 1.7702583261181434080463164711737086e-07),
                      ADDP(-4.8116726436626674386837474600122546e-08),
                      ADDP( 1.8580243803898903932739968462723257e-08),
                      ADDP( 5.9429290158995810274084045450365155e-08),
                      ADDP(-7.3475815263031014338043532409530826e-08),
                      ADDP(-3.6010451791131783412428107502651588e-10),
                      ADDP( 9.1937621332306223709500539337707802e-08),
                      ADDP(-2.3081430789008845105459480604328606e-09),
                      ADDP(-6.3770082710374011644192544857793965e-08),
                      ADDP( 7.0065470710788996937433042385088690e-08),
                      ADDP(-2.6231646223789647162400542135901698e-08),
                      ADDP(-1.6131700179334699793562384581956369e-08)};




    printf("    Surface SHA with a symmetric grid of cells...\n");
    CHARM(err_handler)(err, 1);
    errnum += cmp_arrays(shcs_out->c[0], shcs_out_symm_ref_cnm,
                         ncoeffs, ADDP(1000.0) * CHARM(glob_threshold2));
    errnum += cmp_arrays(shcs_out->s[0], shcs_out_symm_ref_snm,
                         ncoeffs, ADDP(1000.0) * CHARM(glob_threshold2));


    free(f);
    /* ..................................................................... */


    /* Non-symmetric grid */
    /* ..................................................................... */
    f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
    if (f == NULL)
    {
        fprintf(stderr, "Failed to initialize an array to store the "
                        "synthesized signal.\n");
        exit(CHARM_FAILURE);
    }


    /* To get a non-symmetric grid, we can just modify a few latitudes of the
     * previous grid. */
    grd->lat[1] += ADDP(0.005);
    grd->lat[2] += ADDP(0.005);


    /* Do the synthesis in the cells to get a signal to be harmonically
     * analyzed */
    CHARM(shs_cell)(grd, shcs_ref, nmax, f, err);
    CHARM(err_handler)(err, 1);


    CHARM(sha_cell)(grd, f, nmax, CHARM_SHA_CELL_AQ, shcs_out, err);
    CHARM(err_handler)(err, 1);


    REAL shcs_out_nonsymm_ref_cnm[ncoeffs] = {
                        shcs_out->c[0][0], /* We do not want to validate this
                                              coefficient, since it has been
                                              removed from the input
                                              coefficients.  It is non-zero
                                              only due to numerical errors and
                                              is a bit difficult to validate
                                              properly.  By taking the value
                                              from the recovered set of
                                              coefficients, we make sure that
                                              it always passes the comparison
                                              test. */
                        ADDP( 2.5949114551062079705878135906314620e-09),
                        ADDP(-4.7369958669193433952174400891754913e-04),
                        ADDP( 9.3702592218244216098986510570914767e-07),
                        ADDP( 3.8477519901040177469862795750066761e-07),
                        ADDP( 1.1313908167221019692869876259347080e-07),
                        ADDP(-2.9812477077497105376335127695633633e-07),
                        ADDP( 1.5102870713409861750972801221314344e-07),
                        ADDP(-1.7155995858579135021224732348878161e-07),
                        ADDP( 9.5959535591392341593734649322186319e-08),
                        ADDP(-1.8679122280707205971628520814920302e-07),
                        ADDP( 7.9589198894113618629412360233814073e-09),
                        ADDP(-3.8776995469409127955575693829305307e-09),
                        ADDP( 1.9648158371034932352850641701569904e-06),
                        ADDP(-5.0525688635742968987391490079992223e-07),
                        ADDP(-2.8252166278200572889319115744462394e-08),
                        ADDP(-7.2021875351494577914448453238857034e-08),
                        ADDP( 2.5417769308068536958834813669406497e-07),
                        ADDP( 1.2966011402147791870043003692293565e-08),
                        ADDP( 1.2730882267571306299545979116574225e-07),
                        ADDP( 5.1721436116026066462456609725168495e-08),
                        ADDP( 2.3473321670575412895583729427343356e-09),
                        ADDP( 8.8568182827703621618917947653166902e-07),
                        ADDP( 3.3299302321514266526003314517607759e-07),
                        ADDP( 6.1896274276224717733170071573343997e-07),
                        ADDP( 4.9108921236013522608219345296032535e-08),
                        ADDP( 3.0270124264897853665215909038110033e-07),
                        ADDP( 6.6403166304567858010817139961139789e-08),
                        ADDP( 4.3661130067529855897355394997628305e-08),
                        ADDP(-5.9812103507696628239783838855666298e-08),
                        ADDP( 6.9057423262614777722420430102927910e-07),
                        ADDP( 9.3699900901238624476594976514334610e-07),
                        ADDP(-4.0626117327147744490213158103868418e-07),
                        ADDP( 6.6197373380657856121592508174104104e-08),
                        ADDP( 2.0021094904818783994166882241027351e-07),
                        ADDP(-6.3595827736024089380408349363496528e-09),
                        ADDP(-1.1901255983786426292843274917873482e-07),
                        ADDP(-2.7401635543916939674161926523508515e-10),
                        ADDP(-1.7851834034144256123712207928877340e-07),
                        ADDP(-2.7599744360343499089119194126381293e-07),
                        ADDP(-8.6964650037005368456226308211395719e-08),
                        ADDP(-2.4223400779153596100553627392249315e-07),
                        ADDP(-2.0804420824569434553643392084670290e-07),
                        ADDP(-1.9873836205828179776034377587948253e-08),
                        ADDP(-7.4814128405298572795892327030047221e-08),
                        ADDP( 1.5808206386130404328647765801493003e-07),
                        ADDP(-2.3634127260179163274466313487799616e-07),
                        ADDP( 3.3427660433683588216018529715642647e-09),
                        ADDP(-2.9262929724918558880123196387361014e-08),
                        ADDP(-1.2254587194104587940899089676854544e-08),
                        ADDP(-3.9812417954969184805013713832377424e-08),
                        ADDP( 7.2828702659676969873839544001082743e-09),
                        ADDP(-3.0015661085574663021985676209481563e-07),
                        ADDP(-5.4571445941309698612187286611197310e-08),
                        ADDP( 3.9394332783702100779311103660333917e-08),
                        ADDP(-3.0155114268285967173466381853815043e-08),
                        ADDP(-9.7794805669399459753687040895352390e-10),
                        ADDP( 5.3808136744100491899596962536568130e-08),
                        ADDP(-9.0568883816935248588223307604972252e-08),
                        ADDP( 7.8136551377937149646444933979840301e-09),
                        ADDP(-9.5565684380607493459384232249377439e-08),
                        ADDP( 1.4095428252962486243760413436875869e-07),
                        ADDP( 2.6908437074079241849722800819183616e-08),
                        ADDP(-3.4617231927064342658465032010509203e-08),
                        ADDP( 8.7429016526703105762081285702690631e-08),
                        ADDP( 6.7535689277614476977726399377149484e-08)};
    REAL shcs_out_nonsymm_ref_snm[ncoeffs] = {
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 0.0000000000000000000000000000000000e+00),
                        ADDP( 8.7217786425765132008816469082375512e-10),
                        ADDP(-2.2202073185012417043954640958622300e-09),
                        ADDP( 2.4016579556920325100574253355477109e-07),
                        ADDP(-4.4640657073106581599653571502508261e-07),
                        ADDP(-8.1564882353787325656007071874741547e-08),
                        ADDP( 1.5279095390475272446572996641539815e-08),
                        ADDP( 8.0779996845179797205309787615616995e-08),
                        ADDP( 3.8026014093772473361280822837919631e-08),
                        ADDP( 1.8806528769009649537763304619750146e-08),
                        ADDP(-9.3361103074120659700108800245996365e-08),
                        ADDP(-1.3719514939959832324225541159729277e-06),
                        ADDP(-6.0237625952548810353109099665759785e-07),
                        ADDP( 6.1262378364123260749041760757750185e-07),
                        ADDP(-3.0310688651108409544745474097004278e-07),
                        ADDP(-3.2386884338075418960701106972923287e-07),
                        ADDP( 6.5804113783901455501581402336664721e-08),
                        ADDP( 4.4344941015982296881559114811741426e-08),
                        ADDP(-2.9734922519258067466781633583025825e-08),
                        ADDP(-3.8799754616141985453368085337361514e-08),
                        ADDP( 1.3582166651127761111115188206742795e-06),
                        ADDP(-1.9110468344501881351101392037985930e-07),
                        ADDP(-1.9175360210194171263425780826419032e-07),
                        ADDP(-6.5892119998225637446106861760646128e-10),
                        ADDP(-1.8874231032393786451855843778576252e-07),
                        ADDP(-7.8828010082737880568164188476226118e-08),
                        ADDP(-6.5655232583476907114472084801827763e-08),
                        ADDP(-1.1627475673950173728651231091313452e-07),
                        ADDP( 2.8410770318945104791922415727277197e-07),
                        ADDP( 4.3189876281653634926400038336215217e-08),
                        ADDP(-4.1506653386012355563186916401287899e-07),
                        ADDP(-1.0448036203877345211528304936206550e-07),
                        ADDP( 4.1292844375694190934150155548132507e-08),
                        ADDP( 1.0897505963392673960583443324074327e-08),
                        ADDP(-6.0813810069795300538712108467555696e-08),
                        ADDP(-6.0489329500328090702203128741094036e-07),
                        ADDP(-4.7100302378487478352891144610874411e-07),
                        ADDP( 4.7189207715214615699496416422716869e-09),
                        ADDP( 5.8417108346616408294720631329772921e-08),
                        ADDP(-4.4244109428404595826888064693469328e-08),
                        ADDP(-3.7883598924043080535144971578800846e-08),
                        ADDP(-2.0115580634638701909555682725449706e-07),
                        ADDP( 1.3341439107017313830318783577013025e-07),
                        ADDP( 2.4440325563335711833359754966407096e-07),
                        ADDP( 1.7702614350561070322667420553176625e-07),
                        ADDP(-4.8117311904165724048067120866291562e-08),
                        ADDP( 1.8580243309030996392956905314826856e-08),
                        ADDP( 5.9429292052386767655027897929419179e-08),
                        ADDP(-7.3475820464238106524804741792523839e-08),
                        ADDP(-3.6009268228301010394616844515532647e-10),
                        ADDP( 9.1937621277985440347999179428257141e-08),
                        ADDP(-2.3081428593121785156020804482823120e-09),
                        ADDP(-6.3770083346928338050731461845487973e-08),
                        ADDP( 7.0065470712826190589199740057091902e-08),
                        ADDP(-2.6231646232441157762528986778124916e-08),
                        ADDP(-1.6131700179360852659820203099976744e-08)};


    printf("    Surface SHA with a non-symmetric grid of cells...\n");
    errnum += cmp_arrays(shcs_out->c[0], shcs_out_nonsymm_ref_cnm,
                         ncoeffs, ADDP(100.0) * CHARM(glob_threshold2));
    errnum += cmp_arrays(shcs_out->s[0], shcs_out_nonsymm_ref_snm,
                         ncoeffs, ADDP(100.0) * CHARM(glob_threshold2));
    /* ..................................................................... */


    CHARM(shc_free)(shcs_out);
    CHARM(crd_free)(grd);
    free(f);
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs_ref);
    CHARM(err_free)(err);
    free(dda);
    free(dda_ref);


    return errnum;
    /* --------------------------------------------------------------------- */
}
/* ------------------------------------------------------------------------- */

