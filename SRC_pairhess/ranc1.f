       SUBROUTINE RANC1 (IX1, IX2, IX3, IDUM, IFF, R, RANDOM )
       implicit real*8 (a-h,o-z)
!
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
         REAL*8 R(97), RANDOM
         INTEGER IX1, IX2, IX3
         PARAMETER (M1 = 259200, IA1 = 7141, IC1 = 54773, RM1 = 1.D0/259200)
         PARAMETER (M2 = 134456, IA2 = 8121, IC2 = 28411, RM2 = 1.D0/134456)
         PARAMETER (M3 = 243000, IA3 = 4561, IC3 = 51349)
!         IF (IDUM .LT. 0 .OR. IFF .EQ. 0) THEN
         IF( IDUM .LT. 0)THEN
           IFF = 1
           IX1 = MOD (MOD (54773 - IDUM, 259200) * 7141 + 54773, 259200)
           IX2 = MOD (IX1, 134456)
           IX1 = MOD (IX1 * 7141 + 54773, 259200)
           IX3 = MOD (IX1, 243000)
           DO   J=1,97
             IX1 = MOD (IX1 * 7141 + 54773, 259200)
             IX2 = MOD (IX2 * 8121 + 28411, 134456)
             R(J) = (FLOAT (IX1) + FLOAT (IX2) * RM2) * RM1
           END DO
           IDUM = 1
         END IF
         IX1 = MOD (IX1 * 7141 + 54773, 259200)
         IX2 = MOD (IX2 * 8121 + 28411, 134456)
         IX3 = MOD (IX3 * 4561 + 51349, 243000)
         J = (IX3 * 97) / 243000 + 1
         IF (J .GT. 97 .OR. J .LT. 1) THEN
           WRITE (6, *) 'J is out of range, value ', J
           PAUSE
         END IF
         RANDOM = R(J)
         IF (R(J) .GT. 1.D0) RANDOM = R(J) - INT (R(J))
         R(J) = (DBLE (IX1) + DBLE (IX2) * RM2) * RM1
         RETURN
       END
 
