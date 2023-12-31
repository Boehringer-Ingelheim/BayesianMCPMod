# BayesianMCPMod CRAN Comments

## Maintainer

'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'

## Test environments
- aarch64-apple-darwin20, R 4.2.1 (local)
- x86_64-w64-mingw32, (unstable) (2023-12-29 r85751 ucrt)
- r-release-macosx-arm64|4.3.0 (2023-05-18 r84451)
- R-hub Windows Server 2022, R-devel, 64 bit
- R-hub Ubuntu Linux 20.04.1 LTS, R-release
- R-hub Fedora Linux, R-devel

## R CMD check results

### aarch64-apple-darwin20, R 4.2.1 (local)
0 errors √ | 0 warnings √ | 0 notes √

### x86_64-w64-mingw32, (unstable) (2023-12-29 r85751 ucrt)

0 errors | 0 warnings | 1 note

- checking CRAN incoming feasibility
  - Maintainer: 'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'
  - New submission
  - Possibly misspelled words in DESCRIPTION:
    - MCPMod (3:18, 12:23, 14:17, 16:56)
    - Pinheiro (25:25)
    - Schmidli (21:15)
    - al (12:44, 21:27, 25:37)
    - et (12:41, 21:24, 25:34)
    
    => marked spelling is correct

### r-release-macosx-arm64|4.3.0 (2023-05-18 r84451)

0 errors | 0 warnings | 0 notes

### R-hub Windows Server 2022, R-devel, 64 bit

0 errors | 0 warnings | 4 notes

- checking CRAN incoming feasibility
  - Maintainer: 'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'
  - New submission
  - Possibly misspelled words in DESCRIPTION:
    - MCPMod (3:18, 12:23, 14:17, 16:56)
    - Pinheiro (25:25)
    - Schmidli (21:15)
    - al (12:44, 21:27, 25:37)
    - et (12:41, 21:24, 25:34)
    
    => marked spelling is correct
 
- checking HTML version of manual ... NOTE
  - Skipping checking math rendering: package 'V8' unavailable
  
  => R-hub indueced note
 
- checking for non-standard things in the check directory ... NOTE
  - Found the following files/directories:
    - ''NULL''
    
    => R-hub indueced note
 
- checking for detritus in the temp directory ... NOTE
  - Found the following files/directories:
    - 'lastMiKTeXException'
    
    => R-hub indueced note
  
### R-hub Ubuntu Linux 20.04.1 LTS, R-release

0 errors | 0 warnings | 3 notes

- checking CRAN incoming feasibility
  - Maintainer: 'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'
  - New submission
  - Possibly misspelled words in DESCRIPTION:
    - al (12:44, 21:27, 25:37)
    - et (12:41, 21:24, 25:34)
    - MCPMod (3:18, 12:23, 14:17, 16:56)
    - Pinheiro (25:25)
    - Schmidli (21:15)
    
    => marked spelling is correct

- checking examples ... [11s/25s] NOTE
  - Examples with CPU (user + system) or elapsed time > 5s
                          user  system elapsed
    getBootstrapQuantiles 5.78  0.141  13.032
    
    => microbenchmarking the getBootstrapQuantiles() example locally:
         Unit: seconds
                 min         lq       mean     median         uq        max neval
          0.01434709 0.01472304 0.01566722 0.01491748 0.01517197 0.02247079   100
    
- checking HTML version of manual ... NOTE
  - Skipping checking HTML validation: no command 'tidy' found
  - Skipping checking math rendering: package 'V8' unavailable
  
  => R-hub indueced note

### R-hub Fedora Linux, R-devel

0 errors | 0 warnings | 3 notes

- checking CRAN incoming feasibility
  - Maintainer: 'Stephan Wojciekowski <stephan.wojciekowski@boehringer-ingelheim.com>'
  - New submission
  - Possibly misspelled words in DESCRIPTION:
    - MCPMod (3:18, 12:23, 14:17, 16:56)
    - Pinheiro (25:25)
    - Schmidli (21:15)
    - al (12:44, 21:27, 25:37)
    - et (12:41, 21:24, 25:34)
    
    => marked spelling is correct

- checking examples ... [11s/25s] NOTE
  - Examples with CPU (user + system) or elapsed time > 5s
                           user  system  elapsed
    getBootstrapQuantiles 5.601    0.11   12.862
    
    => microbenchmarking the getBootstrapQuantiles() example locally:
         Unit: seconds
                 min         lq       mean     median         uq        max neval
          0.01434709 0.01472304 0.01566722 0.01491748 0.01517197 0.02247079   100
    
- checking HTML version of manual ... NOTE
  - Skipping checking HTML validation: no command 'tidy' found
  - Skipping checking math rendering: package 'V8' unavailable
  
  => R-hub indueced note
     
## From NEWS.md

### BayesianMCPMod 1.0.0 (31-Dec-2023)

- Initial release of the 'BayesianMCPMod' package
- Special thanks to Jana Gierse, Bjoern Bornkamp, Chen Yao, Marius Thoma & Mitchell Thomann for their review and valuable comments
- Thanks to Kevin Kunzmann for R infrastructure support and to Frank Fleischer for methodological support
