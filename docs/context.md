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

- Fantoma geometrico principal de radon/src/phantom.py.
- Parametros iniciales: N = 256 y R = 0.9.
- Mapa de color principal: bone.
- Ruido gaussiano con semilla fija solo si queda documentado.
- Las figuras finales deben ser reproducibles dentro de tesis/ y no depender de ../radon ni ../radon_intr.
