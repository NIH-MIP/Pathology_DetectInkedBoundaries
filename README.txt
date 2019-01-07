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
# (Optional - NEW 1/7/2019) users can elect to manually outline inked borders, 
# which are then automatically registered to un-inked imaging (if provided) and
# annotations automatically saved. This may come in use if the 'auto'
# option fails to accurately detect boundaries
#
# Annotations will be saved to the highest resolution image contained in
# the digital stack (i.e. 40x if full digital image) provided at input.
#
# This program utilizes other freely available toolboxes, see dependencies
#
# General workflow:
#       1. read in image stack, find lowest magnification ratio (smallest
#           sampling of full-res image)
#       2. Option 'auto' ink detection/outlining
#               2a. deconvolve RGB H&E image using Khan et al Classifier
#               2b. Make a mask of tissue sample using H channel
#               2c. Use k-means clustering to detect inked markings from 
#                   remaining channels within tissue sample
#               2d. Grow-shrink morhological operations to determine ROIs
#               2e. Initiate user interface for accept/reject of all 
#                   proposed ROIs
#               2f. Initiate user interface to allow editing of accepted ROIs
#          Option 'manual' ink outlining
#               3a. user prompted for number of ROIs
#               3b. initiate user interface for manual outlining 
#               3c. after each ROI is outlined, user should close figure
#       3. Correlation-based image registration of inked and un-inked
#       4. Write annotation file 
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
#   DetectInkedBoundaries('--Marked','/example/path/to/marked.svs','--Unmarked', '/example/path/to/unmarked.svs','--Method','auto')
#
#   input:
#           1. (required) digital image in CZI or SVS format with inked markings
#           2. (optional) digital image of specimen with removed markings,
#                   if provided, will register low res digital images and register digital
#                   markings to 'clean' specimen
#           3. (optional) method of annotation. default = 'auto'. users who
#                   wish to outline images themselves should use 'manual'
#
## TIPS AND TRICKS #########################################################################################
#
#  - Points can be added and deleted from ROIs after accept/reject stage, this stage is meant to exclude 
#       false positives arising from other inked notes/markings (i.e. "EPE", arrows, etc).
#
#  - If several proposed ROIs are part of the same region, simply add more points to join ROIs as part of
#       Workflow Step 7. A closing operation is performed that will merge all overlapping ROIs.
#
#  - If an ROI is open-ended (i.e. not closed circle), the tool will encompass only ink. To expand to other 
#       regions, all add and/or drag points to encompass the full ROI area
