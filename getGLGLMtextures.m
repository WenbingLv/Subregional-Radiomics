function [textures] = getGLGLMtextures(GLGLM)
% -------------------------------------------------------------------------
% function [textures] = getGLGLMtextures(GLGLM)
% -------------------------------------------------------------------------
% DESCRIPTION:
% This function computes texture features from an input Gray-Level 
% Gap-Length Matrix (GLGLM).
% -------------------------------------------------------------------------
% REFERENCES:
% [1] Galloway, M. M. (1975). Texture analysis using gray level Gap lengths. 
%     Computer Graphics and Image Processing, 4(2), 172–179.
% [2] Chu, A., Sehgal, C. M., & Greenleaf, J. F. (1990). Use of gray value 
%     distribution of Gap lengths for texture analysis. Pattern Recognition
%     Letters, 11(6), 415-419.
% [3] Dasarathy, B. V., & Holder, E. B. (1991). Image characterizations 
%     based on joint gray level-Gap length distributions. Pattern 
%     Recognition Letters, 12(8), 497-502.
% [4] Thibault, G., Fertil, B., Navarro, C., Pereira, S., Cau, P., Levy, 
%     N., Mari, J.-L. (2009). Texture Indexes and Gray Level Size Zone 
%     Matrix. Application to Cell Nuclei Classification. In Pattern 
%     Recognition and Information Processing (PRIP) (pp. 140–145).
% -------------------------------------------------------------------------
% INPUTS:
% - GLGLM: Gray-Level Gap-Length Matrix.
%
% ** 'GLGLM' should be the output from 'getGLGLM.m' **
% -------------------------------------------------------------------------
% OUTPUTS:
% - textures: Struture specifying the values of different GLGLM texture
%             features as defined below.
% -------------------------------------------------------------------------
% AUTHOR(S): Martin Vallieres <mart.vallieres@gmail.com>
% -------------------------------------------------------------------------
% HISTORY:
% - Creation: January 2013
% - Revision: May 2015
% -------------------------------------------------------------------------
% STATEMENT:
% This file is part of <https://github.com/mvallieres/radiomics/>, 
% a package providing MATLAB programming tools for radiomics analysis.
% --> Copyright (C) 2015  Martin Vallieres
%
%    This package is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This package is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this package.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------


% USEFUL MATRICES, VECTORS AND QUANTITIES
sz = size(GLGLM); % Size of GLGLM
nGaps = sum(GLGLM(:));
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLGLM
pg = sum(GLGLM,2)'; % Gray-Level Gap-Number Vector
pr = sum(GLGLM); % Gap-Length Gap-Number Vector


% COMPUTATION OF TEXTURE FEATURES
% 1. Short Gap Emphasis (SGE), Ref.[1]
textures.SGE = (pr*(cVect.^(-2))')/nGaps;

% 2. Long Gap Emphasis (LGE), Ref.[1]
textures.LGE = (pr*(cVect.^2)')/nGaps;

% 3. Gray-Level Fluctuation (GLF), adapted from Ref.[1]
textures.GLF = sum(pg.^2)/nGaps;

% 4. Gap-Length Nonuniformity (GaLN), adapted from Ref.[1]
textures.GaLN = sum(pr.^2)/nGaps;

% 5. Gap Percentage (GP), adapted from Ref.[1]
textures.GP = nGaps/(pr*cVect');

% 6. Low Gray-Level Gap Emphasis (LGGE), Ref.[2]
textures.LGGE = (pg*(rVect.^(-2))')/nGaps;

% 7. High Gray-Level Gap Emphasis (HGGE), Ref.[2]
textures.HGGE = (pg*(rVect.^2)')/nGaps;

% 8. Short Gap Low Gray-Level Emphasis (SGLGE), Ref.[3]
textures.SGLGE = sum(sum(GLGLM.*(rMat.^(-2)).*(cMat.^(-2))))/nGaps;

% 9. Short Gap High Gray-Level Emphasis (SGHGE), Ref.[3]
textures.SGHGE = sum(sum(GLGLM.*(rMat.^2).*(cMat.^(-2))))/nGaps;

% 10. Long Gap Low Gray-Level Emphasis (LGLGE), Ref.[3]
textures.LGLGE = sum(sum(GLGLM.*(rMat.^(-2)).*(cMat.^2)))/nGaps;

% 11. Long Gap High Gray-Level Emphasis (LGHGE), Ref.[3]
textures.LGHGE = sum(sum(GLGLM.*(rMat.^2).*(cMat.^2)))/nGaps;


% New features according to Ref.[4]
GLGLM = GLGLM./nGaps;
pg = sum(GLGLM,2)'; pr = sum(GLGLM);
ug = (pg*rVect')/(sz(1)*sz(2));
ur = (pr*cVect')/(sz(1)*sz(2));

% 12. Gray-Level Variance (GrLV), adapted from Ref.[4]
GrLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        GrLV = GrLV + (GLGLM(g,r)*g-ug)^2;
    end
end
GrLV = GrLV/(sz(1)*sz(2));
textures.GrLV = GrLV;

% 13. Gap-Length Variance (GaLV), adapted from Ref.[4]
GaLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        GaLV = GaLV + (GLGLM(g,r)*r-ur)^2;
    end
end
GaLV = GaLV/(sz(1)*sz(2));
textures.GaLV = GaLV;

end