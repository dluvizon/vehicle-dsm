#include "speed.h"
#include "image.h"

#include <string.h>

#ifdef  __cplusplus
extern "C" void draw_motion_vectors_lum(uint8_t *image, int width, int height,
		int iframe, KLT_FeatureList features);
#else
void draw_motion_vectors_lum(uint8_t *image, int width, int height,
		int iframe, KLT_FeatureList features);
#endif

static const double ref_dist_meters = 1.0;
static const double ref_dist_pixels = 60;

/* reference for 30.15 fps */
static const double frame_interval = 0.033167;

static const double lane_coef[3] = {0.972, 0.933, 0.901};
//static const double lane_coef[3] = {0.915, 0.870, 0.850};

static const int SDEV = 3;
static const double CDEV = 0.5;

static struct hmatrix ipm_mat[3];

static void dump_hmatrix(struct hmatrix *hm)
{
	printf("dump hmatrix:\n");
	printf("%.6f\t%.6f\t%.6f\n", hm->a0, hm->a1, hm->a2);
	printf("%.6f\t%.6f\t%.6f\n", hm->b0, hm->b1, hm->b2);
	printf("%.6f\t%.6f\t%.6f\n", hm->c0, hm->c1, hm->c2);
	printf("\n");
}

int load_ipm_matrix_from_file(const char *filename)
{
	struct hmatrix *hm;
	FILE *file;
	int i;

	file = fopen(filename, "r");
	if (!file) {
		fprintf(stderr, "File not found '%s' ", filename);
		perror("");
		return -1;
	}
	for (i = 0; i < 3; i++) {
		hm = (ipm_mat + i);
		fscanf(file, "%f\t%f\t%f\n", &hm->a0, &hm->a1, &hm->a2);
		fscanf(file, "%f\t%f\t%f\n", &hm->b0, &hm->b1, &hm->b2);
		fscanf(file, "%f\t%f\t%f\n\n", &hm->c0, &hm->c1, &hm->c2);
		dump_hmatrix(hm);
	}

	return 0;
}

void compute_ipm_features(KLT_FeatureList features, int frame_i)
{
	float x;
	float y;
	float w;
	int l;
	int i;

	if ((features->vehicle.lane < 1) ||
			(features->vehicle.lane > GTRUTH_MAX_LANES)) {
		/* This is not expected. */
		printf("Error: features::lane invalid (%d)\n",
				features->vehicle.lane);
		return;
	}
	l = features->vehicle.lane - 1;
	if (l > sizeof(ipm_mat) / sizeof(struct hmatrix)) {
		/* shit happend... */
		printf("Error: lane invalid (%d)\n", l);
		return;
	}

	for (i = 0; i < features->nFeatures; i++) {
		x = features->feature[i]->x;
		y = features->feature[i]->y;
		w = (x * ipm_mat[l].c0) + (y * ipm_mat[l].c1) + ipm_mat[l].c2;
		features->feature[i]->ipm_x[frame_i % SMPSIZE] =
			(x * ipm_mat[l].a0 + y * ipm_mat[l].a1 +
			 ipm_mat[l].a2) / w;
		features->feature[i]->ipm_y[frame_i % SMPSIZE] =
			(x * ipm_mat[l].b0 + y * ipm_mat[l].b1 +
			 ipm_mat[l].b2) / w;
		features->feature[i]->prv_x[frame_i % SMPSIZE] = x;
		features->feature[i]->prv_y[frame_i % SMPSIZE] = y;
	}
}

speed_i compute_velocity_vector(KLT_FeatureList features,
		int frame_i, int frame_j) {

   speed_i speed = {0.0, 0.0, 0.0, 0.0};

   if (frame_i < 1) { return speed; }

   int nfeatures = 0;

   for (int i = 0; i < features->nFeatures; i++)  {

       if (features->feature[i]->val == KLT_TRACKED) {

          double act_x = features->feature[i]->ipm_x[frame_j % SMPSIZE];
          double act_y = features->feature[i]->ipm_y[frame_j % SMPSIZE];
          double prv_x = features->feature[i]->ipm_x[frame_i % SMPSIZE];
          double prv_y = features->feature[i]->ipm_y[frame_i % SMPSIZE];

          if ((act_x == NOT_FILLED) || (act_y == NOT_FILLED) ||
			  (prv_x == NOT_FILLED) || (prv_y == NOT_FILLED)) {
              return speed;
          }
          speed.v_x += (act_x - prv_x);
          speed.v_y += (act_y - prv_y);
          speed.m_x += act_x;
          speed.m_y += act_y;

          nfeatures++;
       }
   }

   speed.v_x = speed.v_x / (double) nfeatures;
   speed.v_y = speed.v_y / (double) nfeatures;
   speed.m_x = speed.m_x / (double) nfeatures;
   speed.m_y = speed.m_y / (double) nfeatures;

   return speed;
}

double metric_calc_velocity_ipm(speed_i speed, unsigned lane)
{
	lane -= 1;
	assert(lane < (sizeof(lane_coef) / sizeof(double)));

	/* Displacement in pixels: */
	double dpixel =
		sqrt((speed.v_x * speed.v_x) + (speed.v_y * speed.v_y));

	/* Distance in meters: */
	double d = (ref_dist_meters * dpixel) / ref_dist_pixels;

	/* Speed in meters per sec: */
	return lane_coef[lane] * d / frame_interval;
}

void outlier_removal(uint8_t *image, int width, int height,
		KLT_FeatureList features, int iframe)
{
   int outliers;
   int nfeatures;
   double sum_x;
   double sum_y;
   uint8_t *tmp1;
   uint8_t *tmp2;
   bool removed;
   int iter = 0;
   int i;

   if (DEBUG_OUTLIER_KLT) {
      tmp1 = (uint8_t *) malloc(width * height * sizeof(uint8_t));
      tmp2 = (uint8_t *) malloc(width * height * sizeof(uint8_t));
   }

   do {
      /* Computing the x-mean and y-mean of all tracked points: */
      iter++;
      nfeatures = 0;
      outliers = 0;

      sum_x = 0.0;
      sum_y = 0.0;
      for (i = 0; i < features->nFeatures; i++)  {
         if (features->feature[i]->val == KLT_TRACKED) {
            double act_x = features->feature[i]->prv_x[(iframe   ) % SMPSIZE];
            double act_y = features->feature[i]->prv_y[(iframe   ) % SMPSIZE];
            double prv_x = features->feature[i]->prv_x[(iframe -1) % SMPSIZE];
            double prv_y = features->feature[i]->prv_y[(iframe -1) % SMPSIZE];

            assert((act_x != NOT_FILLED) && (act_y != NOT_FILLED));
            //if ((prv_x != NOT_FILLED) && (prv_y != NOT_FILLED)) {
            //   goto end;
            //}
            sum_x += (act_x - prv_x);
            sum_y += (act_y - prv_y);
            nfeatures++;
         }
      }

      double mu_x = (sum_x / (double) nfeatures);
      double mu_y = (sum_y / (double) nfeatures);

      /* Computing the x and y standard deviation of all tracked points: */
      sum_x = 0.0;
      sum_y = 0.0;
      for (i = 0; i < features->nFeatures; i++) {
         if (features->feature[i]->val == KLT_TRACKED) {
            double act_x = features->feature[i]->prv_x[(iframe   ) % SMPSIZE];
            double act_y = features->feature[i]->prv_y[(iframe   ) % SMPSIZE];
            double prv_x = features->feature[i]->prv_x[(iframe -1) % SMPSIZE];
            double prv_y = features->feature[i]->prv_y[(iframe -1) % SMPSIZE];
            double d_x = (act_x - prv_x);
            double d_y = (act_y - prv_y);

            sum_x += (d_x - mu_x) * (d_x - mu_x);
            sum_y += (d_y - mu_y) * (d_y - mu_y);
         }
      }

      double dev_x = sqrt(sum_x / (double) (nfeatures - 1));
      double dev_y = sqrt(sum_y / (double) (nfeatures - 1));

      removed = false;

      if (DEBUG_OUTLIER_KLT) {
         memcpy(tmp1, image, width * height * sizeof(uint8_t));
         draw_motion_vectors_lum(tmp1, width, height, iframe, features);
      }

      /* Removing the outliers: */
      for (i = 0; i < features->nFeatures; i++) {
         if (features->feature[i]->val == KLT_TRACKED) {
            double act_x = features->feature[i]->prv_x[(iframe   ) % SMPSIZE];
            double act_y = features->feature[i]->prv_y[(iframe   ) % SMPSIZE];
            double prv_x = features->feature[i]->prv_x[(iframe -1) % SMPSIZE];
            double prv_y = features->feature[i]->prv_y[(iframe -1) % SMPSIZE];
            double d_x = (act_x - prv_x);
            double d_y = (act_y - prv_y);

            if ((fabs(d_x - mu_x) > SDEV * dev_x && (dev_x > CDEV)) ||
                    ((fabs(d_y - mu_y) > SDEV * dev_y) && (dev_y > CDEV))) {

               /* Set the feature as lost: */
               features->feature[i]->val = KLT_OUTLIER;

               /* Increment the number of outliers: */
               outliers++;
	       removed = true;
            }
         }
      }
      if (DEBUG_OUTLIER_KLT && removed) {
         memcpy(tmp2, image, width * height * sizeof(uint8_t));
         draw_motion_vectors_lum(tmp2, width, height, iframe, features);
         write_uchar_to_pgm(tmp1, width, height,
			 "/tmp/slopes/outlier_%05d_%d_tmp1.pgm", iframe, iter);
         write_uchar_to_pgm(tmp2, width, height,
			 "/tmp/slopes/outlier_%05d_%d_tmp2.pgm", iframe, iter);
      }
   } while (outliers > 0);

end:
   if (DEBUG_OUTLIER_KLT) {
      free(tmp1);
      free(tmp2);
   }
}

