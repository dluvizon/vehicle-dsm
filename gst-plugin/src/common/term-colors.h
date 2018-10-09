/** LightWeight Computer Vision Library
 * ------------------------------------
 * @Copyright (C) 2014 SVi Smart Vision
 * @file term-colors.h
 * @brief Defines the UNIX terminal colors.
 * @author Diogo Luvizon <diogo@luvizon.com>
 * @date 16/10/2014
 */

#ifndef __TERM_COLORS_H
#define __TERM_COLORS_H

#define TR_BLK "\e[0;30m"	/**< Text Regular Black. */
#define TR_RED "\e[0;31m"	/**< Text Regular Red. */
#define TR_GRN "\e[0;32m"	/**< Text Regular Green. */
#define TR_YLW "\e[0;33m"	/**< Text Regular Yellow. */
#define TR_BLU "\e[0;34m"	/**< Text Regular Blue. */
#define TR_PUR "\e[0;35m"	/**< Text Regular Purple. */
#define TR_CYN "\e[0;36m"	/**< Text Regular Cyan. */
#define TR_WHT "\e[0;37m"	/**< Text Regular White. */

#define TB_BLK "\e[1;30m"	/**< Text Bold Black. */
#define TB_RED "\e[1;31m"	/**< Text Bold Red. */
#define TB_GRN "\e[1;32m"	/**< Text Bold Green. */
#define TB_YLW "\e[1;33m"	/**< Text Bold Yellow. */
#define TB_BLU "\e[1;34m"	/**< Text Bold Blue. */
#define TB_PUR "\e[1;35m"	/**< Text Bold Purple. */
#define TB_CYN "\e[1;36m"	/**< Text Bold Cyan. */
#define TB_WHT "\e[1;37m"	/**< Text Bold White. */

#define TU_BLK "\e[4;30m"	/**< Text Underline Black. */
#define TU_RED "\e[4;31m"	/**< Text Underline Red. */
#define TU_GRN "\e[4;32m"	/**< Text Underline Green. */
#define TU_YLW "\e[4;33m"	/**< Text Underline Yellow. */
#define TU_BLU "\e[4;34m"	/**< Text Underline Blue. */
#define TU_PUR "\e[4;35m"	/**< Text Underline Purple. */
#define TU_CYN "\e[4;36m"	/**< Text Underline Cyan. */
#define TU_WHT "\e[4;37m"	/**< Text Underline White. */

#define BA_BLK "\e[40m"		/**< Background Black. */
#define BA_RED "\e[41m"		/**< Background Red. */
#define BA_GRN "\e[42m"		/**< Background Green. */
#define BA_YLW "\e[43m"		/**< Background Yellow. */
#define BA_BLU "\e[44m"		/**< Background Blue. */
#define BA_PUR "\e[45m"		/**< Background Purple. */
#define BA_CYN "\e[46m"		/**< Background Cyan. */
#define BA_WHT "\e[47m"		/**< Background White. */

#define TXTRST "\e[0m"		/**< Text Reset. */

#define COLORS_EN

/** Assembly a color text.
 * @param color One of the defined text color.
 * @param text User const text string.
 * @return Assembly color const text string.
 */
#ifdef COLORS_EN
#define cterm(color, text) color text TXTRST
#else
#define cterm(color, text) text
#endif

#endif /* __TERM_COLORS_H */
