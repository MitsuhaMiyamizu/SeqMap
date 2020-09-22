# SeqMap

Use CI to test SeqMap 1.0.13 and mirror the source code stored at the Univ. of Michigan for research purposes only. (http://www-personal.umich.edu/~jianghui/seqmap/)

Ubuntu latest: ![SeqMap](https://github.com/MitsuhaMiyamizu/SeqMap/workflows/seqmap/badge.svg)

macOS Catalina 10.15.6: [![Build Status](https://dev.azure.com/medbioinfo/seqmap/_apis/build/status/macOS%2010.15.6?branchName=master)](https://dev.azure.com/medbioinfo/seqmap/_build/latest?definitionId=1&branchName=master)

Visual Studio 2017: [![Build status](https://ci.appveyor.com/api/projects/status/da4odbki9ikbn05d/branch/master?svg=true)](https://ci.appveyor.com/project/Kiritoalex/seqmap/branch/master)

Visual Studio 2019: [![Build Status](https://dev.azure.com/medbioinfo/seqmap/_apis/build/status/Windows%20Server%202019%20with%20Visual%20Studio%202019?branchName=master)](https://dev.azure.com/medbioinfo/seqmap/_build/latest?definitionId=2&branchName=master)
# Reason why I do this
It seems many research articles only used the old version of SeqMap, albeit there's an updated version of it,
which is 1.0.13, this repository is created not only for the compilation test but also for the binary test.

~~21/09/2020 Update: Use Visual Studio 2017 to compile the code and generate binary which supports x64 on Windows.~~
# Please note:
g++-compiled binary might behaves the way differently compared with the one compiled with Microsoft Visual C++ compiler. 
Therefore it's highly experimental for using MSVC compiled seqmap.
# USE IT AT YOUR OWN RISK, DO NOT USE IT IF YOU HAVEN'T DONE ANY TEST BEFORE.

# References:

Jiang, H., Wong, W.H. (2008) SeqMap: Mapping Massive Amount of Oligonucleotides to the Genome, Bioinformatics, 24(20).

Frankish A, et al (2018) GENCODE reference annotation for the human and mouse genomes.

# To-do
~~fix the warnings prompted during the compilation of MSVC.~~

***Warning***

_THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE._
