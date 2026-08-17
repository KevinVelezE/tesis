# Contexto canonico de la tesis

## Enfoque vigente

La tesis estudia la reconstruccion tomografica bidimensional mediante la transformada de Radon y el metodo de retroproyeccion filtrada (FBP).

Titulo definitivo: Retroproyeccion filtrada en reconstruccion tomografica: un estudio teorico y computacional.

Capitulo 2 definitivo: Transformada de Radon y retroproyeccion filtrada.

Etapa 2: titulo, introduccion y capitulo 2 actualizados conforme a las convenciones canonicas.

Etapa 2.5: figuras pedagogicas vectoriales para Beer-Lambert, geometria de Radon, sinograma y flujo FBP; deben respetar omega como normal, omega_perp como direccion de integracion y rampa |sigma| sin factor 1/2.

Pendiente global: pdfinfo emite `Syntax Warning: Annotation destination array is too short`; no corregido en Etapa 2.

La regularizacion de Tikhonov queda eliminada como linea matematica y computacional de la tesis. No debe desarrollarse una seccion formal, experimentos ni preliminares orientados a Tikhonov.

---

## Convencion de Fourier

Se usa frecuencia en ciclos por unidad, no frecuencia angular.

Para h en una variable:

```text
hhat(sigma) = integral_R h(t) exp(-2 pi i t sigma) dt
h(t)        = integral_R hhat(sigma) exp( 2 pi i t sigma) d sigma
```

Para f en R^2:

```text
fhat(xi) = integral_R2 f(x) exp(-2 pi i x . xi) dx
```

Notacion:

- sigma: frecuencia unidimensional asociada a la variable del detector.
- xi: frecuencia en R^2.
- theta: angulo en [0, pi).

---

## Geometria de Radon

La convencion geometrica canonica es:

```text
omega(theta) = (cos theta, sin theta)
omega_perp(theta) = (-sin theta, cos theta)
```

- omega(theta) es la normal unitaria a la recta.
- omega_perp(theta) es la direccion de integracion sobre la recta.
- t es la distancia firmada de la recta al origen.
- No mezclar esta convencion con integracion sobre todo S^1 sin ajustar factores.

La transformada de Radon se trabaja con theta en [0, pi):

```text
(R f)(t, theta) = integral_R f(t omega(theta) + s omega_perp(theta)) ds
```

---

## Retroproyeccion y FBP

Con theta en [0, pi), se usa directamente:

```text
(R* g)(x) = integral_0^pi g(x . omega(theta), theta) d theta
```

Con la convencion de Fourier anterior, la formula FBP canonica es:

```text
f = R* F_t^{-1}[ |sigma| F_t(R f) ]
```

No se incluye factor 1/2 ni 1/pi en esta normalizacion.

El filtro rampa teorico es:

```text
tau_rampa(sigma) = |sigma|
```

El factor 0.5 presente en codigo previo no es parte de la rampa canonica y no debe incorporarse dentro de tau.

Decision Etapa 3A: la rampa canonica es tau_rampa(sigma) = |sigma|; eliminar scale=0.5 del simbolo. Si aparece una correccion global discreta, debe ser una constante separada, derivada y documentada, nunca oculta dentro de tau.

---

## Simbolos de filtros

Se usa tau(sigma) para el simbolo continuo del filtro, no A(sigma).

Para filtros generales:

```text
g_tau_hat(sigma, theta) = tau(sigma) g_hat(sigma, theta)
f_tau = R* g_tau
```

En la discretizacion se usara M_tau para la matriz diagonal o vector de entradas que representa la multiplicacion por tau en frecuencia.

Decision Etapa 3A: la rampa sera referencia y A-F seran seis filtros adicionales, para siete simbolos en las comparaciones. Todos los filtros se evaluaran sobre la misma realizacion de ruido. RMSE y error relativo se calcularan sobre una region de evaluacion declarada. El caso B con tau_B(0)=0.3 debe discutirse.

El flujo computacional que debe mantenerse consistente es:

```text
FFT -> M_tau -> IFFT -> retroproyeccion
```

---

## Experimentos canonicos previstos

- Fantoma geometrico principal implementado de forma autocontenida en code/tomography.py.
- Parametros iniciales: N = 256 y R = 0.9.
- Mapa de color principal: bone.
- Ruido gaussiano con semilla fija solo si queda documentado.
- Las figuras finales deben ser reproducibles dentro de tesis/ y no depender de ../radon ni ../radon_intr.

---

## Etapa 3B y 3C

Etapa 3B cerrada en commit 43efd429b5c87dab827f98f3cfb7f12de15adc78.

Implementacion numerica reproducible ubicada en tesis/code, con salidas en Figures/Chapter3, Tables/Chapter3 y Results/Chapter3.

Parametros canonicos del capitulo 3:

```text
N = 256
R = 0.9
N_t = 256
N_fft = 1024
M = 90
Delta theta = 2 grados
seed = 42
noise_level = 0.05
circle = True
```

Decision Etapa 3C: titulo del capitulo 3 cambiado a Implementacion numerica y analisis de la retroproyeccion filtrada.

El filtrado discreto debe explicarse con zero-padding centrado P, recorte C, DFT F_Nfft y multiplicador M_tau evaluado sobre fftfreq. N_fft es longitud auxiliar y no agrega datos ni resolucion.

La tabla de metricas de Chapter3 incluye RMSE, error relativo en norma discreta ell^2 y razon de medias. RMSE y error relativo son metricas principales; la razon de medias es solo diagnostico de amplitud y no factor de escala.

Se agrega referencia formal a scikit-image solamente para justificar el uso de skimage.transform.radon en la generacion del sinograma.

---

## Etapa 4A bibliografica

Etapa 4A incorpora fortalecimiento bibliografico sin cambiar formulas, estructura, codigo, figuras, tablas, metricas ni resultados.

Fuentes nuevas citables agregadas:

- helgason1999 para geometria integral, transformada de Radon, teorema de corte e inversion.
- beckmann_iske2017 para estimaciones de error en FBP.
- beckmann_nickel2025 para diseno contemporaneo de funciones filtro optimizadas en FBP.

Restricciones bibliograficas vigentes:

- Los filtros A-F de la tesis son exploratorios y no reproducen los filtros optimizados de Beckmann-Nickel.
- Beylkin se cita solo para distinguir formulaciones formales de transformada de Radon discreta; el codigo no implementa exactamente esa transformada.
- barrett_myers2013 queda sin citar mientras no se fije inequívocamente la edicion usada.
- No citar Davison, Stein-Wainger ni Shepp-Logan salvo que se agregue contenido explicito que lo requiera.
- Quinto y Frikel-Quinto quedan reservados para conclusiones futuras solo si se propone angulo limitado o datos incompletos.

---

## Etapa 4B cierre editorial

Etapa 4B incorpora resumen real, conclusiones, bibliografia en espanol y encabezado abreviado del capitulo 3 sin cambiar formulas, codigo, figuras, tablas, metricas ni resultados.

Decisiones editoriales:

- No se agrega abstract en ingles ni palabras clave.
- Las conclusiones se incluyen antes de la bibliografia.
- La bibliografia debe mostrarse como Bibliografia mediante titulo explicito de biblatex.
- El titulo visible y el indice del capitulo 3 conservan el nombre completo; solo la marca de encabezado se abrevia como Implementacion numerica y analisis de FBP.
- Quinto1993 y frikel_quinto2013 se citan solo en trabajo futuro por datos incompletos y angulo limitado.

Interpretacion de cierre que debe preservarse:

- La rampa canonica es tau_rampa(sigma)=|sigma|, sin factores adicionales.
- La rampa obtiene el menor error sin ruido; E y D ocupan las posiciones siguientes con resultados cercanos entre si.
- C obtiene el menor error con ruido solo bajo la configuracion y realizacion ruidosa estudiadas.
- No se afirma superioridad universal de ningun filtro.
- Las limitaciones son fantoma sintetico, geometria paralela bidimensional, una realizacion de ruido, filtros exploratorios y ausencia de datos clinicos o tridimensionales.

---

## Etapa 5A

- La malla espacial es semiabierta en [-1,1), con h=2/N.
- El origen espacial y el origen del detector corresponden al indice 128.
- CFG.circle se conecta directamente con el parametro circle de skimage.transform.radon.
- El hash reproducible de metricas es 5c898498df56fe763465ca6729afd7eb60a9e9d3613f45421032ff357c8c2802.
- Sin ruido, la rampa ocupa la primera posicion, seguida por E y D.
- Con ruido, C ocupa la primera posicion solo para la configuracion y realizacion estudiadas.
