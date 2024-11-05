
function Img = acode_main_retin_seg(im, bw_mask)

addpath(genpath('./functions'));

bw_mask = logical(bw_mask);

im = mat2gray(im).*mat2gray(bw_mask);
im = imcomplement(im);     
im = im2double(im);
%%
DEG_NUM = 12;
LEN_c = 11;
LEN_o = 11;
LEN_diff = 7;
%
ic1 = reconstruction_by_dilation(im,LEN_c,DEG_NUM); 
io1 = min_openings(im,LEN_o,DEG_NUM);
iv = mat2gray(ic1-io1);
imDiff = smooth_cross_section(iv,LEN_diff,LEN_c);
imL = reconstruction_by_dilation(imDiff,LEN_c,DEG_NUM);
imF = reconstruction_by_erosion(imL,LEN_c,DEG_NUM);

Img = imF;

end








