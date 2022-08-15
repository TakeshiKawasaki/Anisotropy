This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].<br>
[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg

## Obtaining the codes via Git (available in our clusters) <br>
 `git clone https://github.com/TakeshiKawasaki/Anisotropy` repo <br>
 For update:<br>
 `git pull` 
 
## Model <br>
$$
\begin{aligned}
&\sigma(\hat{\mathbf{u}}_{1}, \hat{\mathbf{u}}_{2}, \hat{\mathbf{r}}) \\
&\quad=\sigma_{0}(1-\frac{1}{2} \chi\{\frac{(\hat{\mathbf{r}} \cdot \hat{\mathbf{u}}_{1}+\hat{\mathbf{r}} \cdot \hat{\mathbf{u}}_{2})^{2}}{1+\chi(\mathbf{u}_{1} \cdot u_{2})}+\frac{(\hat{\mathbf{r}} \cdot \hat{\mathbf{u}}_{1}-\hat{\mathbf{r}} \cdot \hat{\mathbf{u}}_{2})^{2}}{1-\chi(\mathbf{u}_{1} \cdot \hat{\mathbf{u}}_{2})}\})^{-1 / 2}
\end{aligned}
$$
