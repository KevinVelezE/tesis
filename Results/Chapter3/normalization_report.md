# Auditoría de normalización

Rampa implementada: `tau_rampa(sigma)=|sigma|`; no se usa `scale=0.5`.

El sinograma se calcula con `skimage.transform.radon(..., circle=True)`.

N_t = 256

N_fft = 1024

Delta x = 0.0078125

Delta t = 0.0078125

Delta sigma = 0.125

Delta theta = 0.0349065850399 rad (2 grados)

El sinograma de longitud `N_t` se rellena de forma centrada hasta `N_fft`, se transforma con `fft(..., axis=0)`, se multiplica como `tau[:, None] * G`, se retorna con `ifft(..., axis=0)` y se recorta al intervalo detector original.

La grilla de frecuencias es `scipy.fft.fftfreq(N_fft, d=Delta t)`. Los símbolos se multiplican en el orden nativo de `fftfreq`; `fftshift` se usa solo para graficar.

La retroproyección evalúa `t=x cos(theta)-y sin(theta)` por interpolación unidimensional en detector y aplica cuadratura angular con `Delta theta`. El signo menos alinea `y` cartesiano con las filas de imagen usadas por `skimage.transform.radon`; la prueba de orientación con fantoma asimétrico verifica esta elección.

Ruido: `sigma_noise = 0.05 * max(abs(g)) = 0.170757766058`; la misma realización ruidosa se usa para todos los filtros.

Valores en frecuencia cero: `{"A": 0.0, "B": 0.3, "C": 0.0, "D": 0.0, "E": 0.0, "F": 0.0, "ramp": 0.0}`. En particular, `tau_B(0)=0.3`.

Centros de la prueba de orientación: `{"phantom_center": [-0.0792237755847951, 0.11746809515107226], "ramp_reconstruction_center": [-0.05857579294066695, 0.09802551214006061]}`. Los signos coinciden, sin rotación/reflexión detectada.

No se aplica corrección empírica ni normalización por reconstrucción. La razón media de la rampa sin ruido en el disco de evaluación es 0.976683; se conserva solo como diagnóstico, no como factor de escala.

Hash reproducible de métricas redondeadas: `5c898498df56fe763465ca6729afd7eb60a9e9d3613f45421032ff357c8c2802`.
