#include "pgrp_dat.h"

int  lat_C2_sgrp[NC2_sgrp]={
   MONOCLINIC_P, /*  1  */
   MONOCLINIC_P, /*  2  */
   MONOCLINIC_A  /*  3  */
};
int  lat_Cs_sgrp[NCs_sgrp]={
   MONOCLINIC_P, /*  1  */
   MONOCLINIC_P, /*  2  */
   MONOCLINIC_A, /*  3  */
   MONOCLINIC_A  /*  4  */
};
int  lat_C2h_sgrp[NC2h_sgrp]={
   MONOCLINIC_P, /*  1  */
   MONOCLINIC_P, /*  2  */
   MONOCLINIC_A, /*  3  */
   MONOCLINIC_P, /*  4  */
   MONOCLINIC_P, /*  5  */
   MONOCLINIC_A  /*  6  */
};

char  *comnt_C2_sgrp[NC2_sgrp]={
  "3 (P 2) [unique axis c]" ,
  "4 (P 21) [unique axis c]",
  "5 (C 2) [unique axis c] cell choice 1"
};
char  *comnt_Cs_sgrp[NCs_sgrp]={
  "6 (P m) [unique axis c]",
  "7 (P c) [unique axis c] cell choice 1",
  "8 (C m) [unique axis c] cell choice 1",
  "9 (C c) [unique axis c] cell choice 1"
};
char  *comnt_C2h_sgrp[NC2h_sgrp]={
  "10 (P 2/m) [unique axis c]"  ,
  "11 (P 21/m) [unique axis c]" ,
  "12 (C 2/m) [unique axis c] cell choice 1" ,
  "13 (P 2/c) [unique axis c] cell choice 1"  ,
  "14 (P 21/c) [unique axis c] cell choice 1" ,
  "15 (C 2/c) [unique axis c] cell choice 1"
};

double rC2_sgrp[NC2_sgrp][2*3]={
  { /*  1  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  },
  { /*  2  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.5000  /*    2  */
  },
  { /*  3  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  }
};
double rCs_sgrp[NCs_sgrp][2*3]={
  { /*  1  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  },
  { /*  2  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.0000  /*    2  */
  },
  { /*  3  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  },
  { /*  4  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.0000  /*    2  */
  }
};
double rC2h_sgrp[NC2h_sgrp][4*3]={
  { /*  1  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.0000,  0.0000,  0.0000  /*    4  */
  },
  { /*  2  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.5000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.0000,  0.0000,  0.5000  /*    4  */
  },
  { /*  3  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.0000,  0.0000,  0.0000  /*    4  */
  },
  { /*  4  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.5000,  0.0000,  0.0000  /*    4  */
  },
  { /*  5  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.5000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.5000,  0.0000,  0.5000  /*    4  */
  },
  { /*  6  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.5000,  0.0000,  0.0000  /*    4  */
  }
};

/* the same as before only B-centered for monoclinic */

int  lat_B_C2_sgrp[NC2_sgrp]={
   MONOCLINIC_P, /*  1  */
   MONOCLINIC_P, /*  2  */
   MONOCLINIC_B  /*  3  */
};
int  lat_B_Cs_sgrp[NCs_sgrp]={
   MONOCLINIC_P, /*  1  */
   MONOCLINIC_P, /*  2  */
   MONOCLINIC_B, /*  3  */
   MONOCLINIC_B  /*  4  */
};
int  lat_B_C2h_sgrp[NC2h_sgrp]={
   MONOCLINIC_P, /*  1  */
   MONOCLINIC_P, /*  2  */
   MONOCLINIC_B, /*  3  */
   MONOCLINIC_P, /*  4  */
   MONOCLINIC_P, /*  5  */
   MONOCLINIC_B  /*  6  */
};

char  *comnt_B_C2_sgrp[NC2_sgrp]={
  "3 (P 2) [unique axis c]" ,
  "4 (P 21) [unique axis c]",
  "5 (C 2) [unique axis c] cell choice 2"
};
char  *comnt_B_Cs_sgrp[NCs_sgrp]={
  "6 (P m) [unique axis c]",
  "7 (P c) [unique axis c] cell choice 1",
  "8 (C m) [unique axis c] cell choice 2",
  "9 (C c) [unique axis c] cell choice 2"
};
char  *comnt_B_C2h_sgrp[NC2h_sgrp]={
  "10 (P 2/m) [unique axis c]"  ,
  "11 (P 21/m) [unique axis c]" ,
  "12 (C 2/m) [unique axis c] cell choice 2" ,
  "13 (P 2/c) [unique axis c] cell choice 1"  ,
  "14 (P 21/c) [unique axis c] cell choice 1" ,
  "15 (C 2/c) [unique axis c] cell choice 2"
};

double rB_C2_sgrp[NC2_sgrp][2*3]={
  { /*  1  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  },
  { /*  2  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.5000  /*    2  */
  },
  { /*  3  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  }
};
double rB_Cs_sgrp[NCs_sgrp][2*3]={
  { /*  1  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  },
  { /*  2  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.0000  /*    2  */
  },
  { /*  3  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000  /*    2  */
  },
  { /*  4  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.5000,  0.0000  /*    2  */
  }
};
double rB_C2h_sgrp[NC2h_sgrp][4*3]={
  { /*  1  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.0000,  0.0000,  0.0000  /*    4  */
  },
  { /*  2  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.5000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.0000,  0.0000,  0.5000  /*    4  */
  },
  { /*  3  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.0000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.0000,  0.0000,  0.0000  /*    4  */
  },
  { /*  4  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.5000,  0.0000,  0.0000  /*    4  */
  },
  { /*  5  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.0000,  0.5000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.5000,  0.0000,  0.5000  /*    4  */
  },
  { /*  6  */
     0.0000,  0.0000,  0.0000, /*    1  */
     0.5000,  0.5000,  0.0000, /*    2  */
     0.0000,  0.0000,  0.0000, /*    3  */
     0.5000,  0.5000,  0.0000  /*    4  */
  }
};
