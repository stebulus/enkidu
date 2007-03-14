\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{enkidudoc}[2007/03/14]
\LoadClass{article}
\RequirePackage{euler}
\RequirePackage{palatino}
\RequirePackage{graphics}
\RequirePackage{hyperref}
\RequirePackage[tiny]{titlesec}
\titleformat*{\subsection}{\rmfamily}
\RequirePackage{tocloft}
\renewcommand{\cfttoctitlefont}{\hfill\rm\bf}
\renewcommand{\cftaftertoctitle}{\hfill}
\renewcommand{\cftsecfont}{\rm}
\renewcommand{\cftsecpagefont}{\rm}
\setlength{\cftbeforetoctitleskip}{0pt}
\setlength{\cftbeforesecskip}{0pt}
\RequirePackage{fancyhdr}
\pagestyle{fancy}
\lhead{} \chead{} \rhead{} \renewcommand{\headrulewidth}{0pt}
\newcommand{\websiteurl}[1]{\def\@websiteurl{#1}}
\newcommand{\websitetext}[1]{\def\@websitetext{#1}}
\newcommand{\websitehref}{\href{\@websiteurl}{\texttt{\@websitetext}}}
\lfoot{\footnotesize \@author\ $\cdot$ \@date\ $\cdot$ \websitehref}
    \cfoot{} \rfoot{\footnotesize\thepage}
\renewcommand{\footrulewidth}{0pt}
\renewcommand{\maketitle}{\centerline{\textbf{\@title}}%
    \vspace{\baselineskip}\hrule\vspace{\baselineskip}}
\newcommand{\fig}[1]{\input{#1-fig.tex}}
\newcommand{\src}[1]{\input{#1-src.tex}}
\newcommand{\figsrc}[2]{%
    \rightline{\raisebox{#2}[0pt][0pt]{\fig{#1}}}%
    \vspace{-\baselineskip}% back up over the empty line
    \src{#1}}
