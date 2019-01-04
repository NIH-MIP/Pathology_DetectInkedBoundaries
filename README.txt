## DetectInkedBoundaries ####################################################################################
#
# Stephanie Harmon
# stephanie.harmon@nih.gov
# Leidos Biomedical Research
# Molecular Imaging Branch, National Cancer Institute
# National Institutes of Health
# December 2018
# 
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
## DESCRIPTION ##############################################################################################
#
# This program takes in digital pathology (.czi/.svs) files with inked
# markings, automatically detects boundaries of ink, provides user
# interface for editing of detected boundaries, and writes to annotation
# file (.cz/.xml). 
#
# (Optional) Additionally, users can provide corresponding un-inked
# digital imaging, which will be automatically registered to inked imaging
# and annotations saved
#
# Annotations will be saved to the highest resolution image contained in
# the digital stack (i.e. 40x if full digital image) provided at input.
#
# This program utilizes other freely available toolboxes, see dependencies
#
# General workflow:
#       1. read in image stack, find lowest magnification ratio (smallest
#           sampling of full-res image)
#       2. deconvolve RGB H&E image using Khan et al Random Forest Classifier
#       3. Make a mask of tissue sample using H channel
#       4. Use k-means clustering to detect inked markings from remaining
#           channels within tissue sample
#       5. Grow-shrink morhological operations to determine ROIs
#       6. Initiate user interface for accept/reject of all proposed ROIs
#       7. Initiate user interface to allow editing of accepted ROIs
#       8. Correlation-based image registration of inked and un-inked
#       9. Write annotation file 
#           (.xml if supplied .svs or .cz if supplied .czi)
#
## DEPENDENCIES #############################################################################################
#
# - local installation of MATLAB 2018b
#
# - MATLAB Image Processing Toolbox
#
# - Bio-Formats toolbox for MATLAB
#     https://docs.openmicroscopy.org/bio-formats/5.3.4/users/matlab/index.html
#
# - color decon toolbox from Khan et al 
#     https://github.com/lun5/color-deconvolution/tree/master/stain_normalisation_toolbox
#
# - polygon decimate function from 
#     https://www.mathworks.com/matlabcentral/fileexchange/34639-decimate-polygon   
#
## USAGE ###################################################################################################
#   
#   DetectInkedBoundaries('--Marked','/example/path/to/marked.svs','--Unmarked', '/example/path/to/unmarked.svs')
#
#   input:
#           1. (required) digital image in CZI or SVS format with inked markings
#           2. (optional) digital image of specimen with removed markings,
#                   if provided, will register low res digital images and register digital
#                   markings to 'clean' specimen 
#
#
