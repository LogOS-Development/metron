//! Tests for metron.

#![cfg_attr(test, allow(clippy::module_inception))]

#[cfg(test)]
mod tests {
    use crate::dim;
    use crate::prefix::SiPrefix;
    use crate::si::*;
    use crate::{
        Acceleration, Area, Capacitance, Charge, Density, Dimensionless, Energy, Force, Frequency,
        GravitationalParameter, Inductance, Kelvins, Kilograms, Kilometers, MagneticFlux,
        MagneticFluxDensity, MassFlowRate, Meters, Millimeters, Milliseconds, MomentOfInertia,
        Nanoteslas, Power, Pressure, Resistance, Seconds, Velocity, Voltage, Volume, Wavenumber,
    };
    use crate::{ConvertPrefix, Quantity, TensorQuantity, UnitName, VectorQuantity};
    use approx::assert_relative_eq;
    use core::marker::PhantomData;
    use core::ops::{Div, Mul};
    use num_traits::Zero;

    // --- Scalar tests ---

    #[test]
    fn scalar_construction() {
        let d = Meters::new(10.0);
        assert_relative_eq!(d.value, 10.0);
        assert_relative_eq!(*d.value(), 10.0);
        assert_relative_eq!(d.into_value(), 10.0);
    }

    #[test]
    fn scalar_add_sub() {
        assert_relative_eq!((Meters::new(3.0) + Meters::new(7.0)).value, 10.0);
        assert_relative_eq!((Meters::new(7.0) - Meters::new(3.0)).value, 4.0);
    }

    #[test]
    fn scalar_mul_div_derives_units() {
        let vel: Velocity = Meters::new(10.0) / Seconds::new(2.0);
        assert_relative_eq!(vel.value, 5.0);
        assert_eq!(dim::Velocity::NAME, "m/s");
    }

    #[test]
    fn scalar_times_raw() {
        assert_relative_eq!((Meters::new(3.0) * 2.0).value, 6.0);
    }

    #[test]
    fn acceleration_chain() {
        let vel: Velocity = Meters::new(10.0) / Seconds::new(2.0);
        let acc: Acceleration = vel / Seconds::new(5.0);
        assert_relative_eq!(acc.value, 1.0);
    }

    #[test]
    fn force_from_mass_x_accel() {
        let mass = Kilograms::new(2.0);
        let acc: Acceleration = Meters::new(10.0) / (Seconds::new(1.0) * Seconds::new(1.0));
        let f: Force = mass * acc;
        assert_relative_eq!(f.value, 20.0);
        assert_eq!(dim::Force::NAME, "N");
    }

    #[test]
    fn prefix_conversion() {
        let km = Kilometers::new(1000.0);
        assert_relative_eq!(km.value, 1000.0);
        assert_relative_eq!(km.in_prefix(SiPrefix::Kilo), 1.0);
    }

    #[test]
    fn display() {
        assert_eq!(format!("{}", Meters::new(42.0)), "42 m");
        let vel: Velocity = Meters::new(10.0) / Seconds::new(2.0);
        assert_eq!(format!("{vel}"), "5 m/s");
    }

    #[test]
    fn default_zero() {
        assert!(Meters::<f64>::default().value == 0.0);
        assert!(Velocity::<f64>::default().value == 0.0);
    }

    #[test]
    fn scalar_neg() {
        assert_relative_eq!((-Meters::new(5.0)).value, -5.0);
    }

    #[test]
    fn scalar_sub() {
        assert_relative_eq!((Seconds::new(10.0) - Seconds::new(3.0)).value, 7.0);
    }

    #[test]
    fn scalar_add_assign() {
        let mut val = Meters::new(5.0);
        val += Meters::new(3.0);
        assert_relative_eq!(val.value, 8.0);
    }

    #[test]
    fn scalar_sub_assign() {
        let mut val = Meters::new(10.0);
        val -= Meters::new(4.0);
        assert_relative_eq!(val.value, 6.0);
    }

    #[test]
    fn scalar_div_raw() {
        assert_relative_eq!((Meters::new(10.0) / 2.0).value, 5.0);
    }

    #[test]
    fn scalar_mul_assign_raw() {
        let mut val = Meters::new(3.0);
        val *= 2.0;
        assert_relative_eq!(val.value, 6.0);
    }

    #[test]
    fn scalar_map_preserves_unit() {
        let dist = Meters::new(-5.0);
        let abs_d = dist.map(f64::abs);
        assert_relative_eq!(abs_d.value, 5.0);
        assert_eq!(dim::Meter::NAME, "m");
    }

    #[test]
    fn scalar_display_all_units() {
        assert_eq!(format!("{}", Force::new(10.0)), "10 N");
        assert_eq!(format!("{}", Energy::new(5.0)), "5 J");
        assert_eq!(format!("{}", Power::new(100.0)), "100 W");
        assert_eq!(format!("{}", Pressure::new(101325.0)), "101325 Pa");
        assert_eq!(format!("{}", Area::new(12.0)), "12 m²");
        assert_eq!(format!("{}", Volume::new(3.0)), "3 m³");
        assert_eq!(format!("{}", Density::new(1000.0)), "1000 kg/m³");
        assert_eq!(format!("{}", Frequency::new(60.0)), "60 Hz");
        let charge_str = format!("{}", Charge::new(1.6e-19));
        assert!(charge_str.ends_with(" C"));
        assert_eq!(format!("{}", Voltage::new(120.0)), "120 V");
        assert_eq!(format!("{}", Resistance::new(50.0)), "50 Ω");
        let cap_str = format!("{}", Capacitance::new(1e-6));
        assert!(cap_str.ends_with(" F"));
        assert_eq!(format!("{}", Inductance::new(0.1)), "0.1 H");
        assert_eq!(format!("{}", MagneticFlux::new(0.5)), "0.5 Wb");
        assert_eq!(format!("{}", MagneticFluxDensity::new(0.05)), "0.05 T");
        assert_eq!(format!("{}", MomentOfInertia::new(100.0)), "100 kg·m²");
        assert_eq!(format!("{}", Wavenumber::new(5.0)), "5 1/m");
        assert_eq!(format!("{}", MassFlowRate::new(2.5)), "2.5 kg/s");
        assert_eq!(format!("{}", Dimensionless::new(0.5)), "0.5 ");
        assert_eq!(format!("{}", Kelvins::new(300.0)), "300 K");
    }

    #[test]
    fn scalar_zero_impl() {
        let z = Meters::<f64>::zero();
        assert!(z.is_zero());
        assert_relative_eq!(z.value, 0.0);
        let nz = Meters::new(1.0);
        assert!(!nz.is_zero());
    }

    #[test]
    fn scalar_deref() {
        let d = Meters::new(10.0);
        assert_relative_eq!(*d, 10.0);
    }

    #[test]
    fn scalar_deref_mut() {
        let mut val = Meters::new(5.0);
        *val += 3.0;
        assert_relative_eq!(val.value, 8.0);
    }

    #[test]
    fn convert_prefix() {
        let meters = Meters::new(1000.0);
        assert_relative_eq!(meters.in_prefix(SiPrefix::Kilo), 1.0);
        assert_relative_eq!(meters.in_prefix(SiPrefix::None), 1000.0);
        let same = meters.convert_to(SiPrefix::Kilo);
        assert_relative_eq!(same.value, 1000.0);
    }

    #[test]
    fn si_prefix_display() {
        assert_eq!(format!("{}", SiPrefix::Kilo), "k");
        assert_eq!(format!("{}", SiPrefix::Milli), "m");
        assert_eq!(format!("{}", SiPrefix::None), "");
        assert_eq!(format!("{}", SiPrefix::Micro), "µ");
        assert_eq!(format!("{}", SiPrefix::Mega), "M");
    }

    #[test]
    fn si_prefix_scale() {
        assert_relative_eq!(SiPrefix::Kilo.scale(), 1.0e3);
        assert_relative_eq!(SiPrefix::Milli.scale(), 1.0e-3);
        assert_relative_eq!(SiPrefix::None.scale(), 1.0);
        assert_relative_eq!(SiPrefix::Yocto.scale(), 1.0e-24);
        assert_relative_eq!(SiPrefix::Yotta.scale(), 1.0e24);
    }

    #[test]
    fn unit_mul_type_level() {
        let _: <dim::Meter as Mul<dim::Second>>::Output = crate::Unit(PhantomData);
    }

    #[test]
    fn unit_div_type_level() {
        let _: <dim::Meter as Div<dim::Second>>::Output = crate::Unit(PhantomData);
    }

    #[test]
    fn gm_type() {
        let mu: GravitationalParameter = Quantity::new(3.986004415e14);
        assert_relative_eq!(mu.value, 3.986004415e14);
        assert_eq!(dim::GravitationalParameter::NAME, "m³/s²");
    }

    // --- Vector tests ---

    #[test]
    fn vector_norm() {
        let pos = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(3.0, 4.0, 0.0));
        assert_relative_eq!(pos.norm().value, 5.0);
    }

    #[test]
    fn vector_normalize() {
        let pos = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(3.0, 4.0, 0.0));
        let dir = pos.normalize();
        assert_relative_eq!(dir.vector.x, 0.6);
        assert_relative_eq!(dir.vector.y, 0.8);
    }

    #[test]
    fn vector_add() {
        let v1 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        let v2 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        let c = v1 + v2;
        assert_relative_eq!(c.vector.x, 1.0);
        assert_relative_eq!(c.vector.y, 1.0);
    }

    #[test]
    fn vector_cross() {
        let xaxis =
            VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        let yaxis =
            VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(0.0, 1.0, 0.0));
        let zaxis = xaxis.cross(&yaxis);
        assert_relative_eq!(zaxis.vector.z, 1.0);
    }

    #[test]
    fn scalar_times_vector() {
        let scale = Meters::new(2.0);
        let dir = VectorQuantity::<f64, 3, dim::Dimensionless>::new(nalgebra::Vector3::new(
            1.0, 2.0, 3.0,
        ));
        let p = scale * dir;
        assert_relative_eq!(p.vector.x, 2.0);
    }

    #[test]
    fn vector_norm_squared() {
        let pos = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(3.0, 4.0, 0.0));
        assert_relative_eq!(pos.norm_squared().value, 25.0);
    }

    #[test]
    fn vector_dot_same_unit() {
        let v1 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 2.0, 3.0));
        let v2 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(4.0, 5.0, 6.0));
        let d = v1.dot(&v2);
        assert_relative_eq!(d.value, 32.0);
    }

    #[test]
    fn vector_dot_different_unit() {
        let pos = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 0.0, 0.0));
        let vel =
            VectorQuantity::<f64, 3, dim::Velocity>::new(nalgebra::Vector3::new(10.0, 0.0, 0.0));
        let d = pos.dot(&vel);
        assert_relative_eq!(d.value, 10.0);
    }

    #[test]
    fn vector_sub() {
        let v1 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(5.0, 3.0, 1.0));
        let v2 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(2.0, 1.0, 0.0));
        let c = v1 - v2;
        assert_relative_eq!(c.vector.x, 3.0);
        assert_relative_eq!(c.vector.y, 2.0);
        assert_relative_eq!(c.vector.z, 1.0);
    }

    #[test]
    fn vector_neg() {
        let v1 = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, -2.0, 3.0));
        let neg = -v1;
        assert_relative_eq!(neg.vector.x, -1.0);
        assert_relative_eq!(neg.vector.y, 2.0);
        assert_relative_eq!(neg.vector.z, -3.0);
    }

    #[test]
    fn vector_div_raw() {
        let v1 =
            VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(10.0, 20.0, 30.0));
        let result = v1 / 2.0;
        assert_relative_eq!(result.vector.x, 5.0);
        assert_relative_eq!(result.vector.y, 10.0);
        assert_relative_eq!(result.vector.z, 15.0);
    }

    #[test]
    fn vector_display() {
        let vq = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 2.0, 3.0));
        let str_v = format!("{vq}");
        assert!(str_v.contains("m"));
    }

    #[test]
    fn vector_into_vector() {
        let vq = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 2.0, 3.0));
        let raw = vq.into_vector();
        assert_relative_eq!(raw.x, 1.0);
    }

    #[test]
    fn vector_accessors() {
        let vq = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(7.0, 8.0, 9.0));
        assert_relative_eq!(vq.vector().x, 7.0);
        assert_relative_eq!(vq.value().x, 7.0);
        assert_relative_eq!(vq.raw().x, 7.0);
    }

    #[test]
    fn vector_default_zero() {
        let vq = VectorQuantity::<f64, 3, dim::Meter>::default();
        assert_relative_eq!(vq.vector.norm(), 0.0);
    }

    // --- Tensor tests ---

    #[test]
    fn tensor_identity() {
        let t: TensorQuantity<f64, 3, 3, dim::MomentOfInertia> = TensorQuantity::identity();
        assert_relative_eq!(t.matrix[(0, 0)], 1.0);
    }

    #[test]
    fn tensor_times_vector() {
        let mat: TensorQuantity<f64, 3, 3, dim::Dimensionless> = TensorQuantity::identity();
        let vq = VectorQuantity::<f64, 3, dim::Meter>::new(nalgebra::Vector3::new(1.0, 2.0, 3.0));
        let r = mat * vq;
        assert_relative_eq!(r.vector.x, 1.0);
    }

    #[test]
    fn tensor_transpose() {
        let mat = TensorQuantity::<f64, 2, 3, dim::Dimensionless>::new(
            nalgebra::SMatrix::from_element(1.0),
        );
        let t = mat.transpose();
        assert_relative_eq!(t.matrix[(0, 0)], 1.0);
    }

    #[test]
    fn tensor_sub() {
        let t1: TensorQuantity<f64, 2, 2, dim::Dimensionless> =
            TensorQuantity::new(nalgebra::SMatrix::from_element(5.0));
        let t2: TensorQuantity<f64, 2, 2, dim::Dimensionless> =
            TensorQuantity::new(nalgebra::SMatrix::from_element(3.0));
        let c = t1 - t2;
        assert_relative_eq!(c.matrix[(0, 0)], 2.0);
    }

    #[test]
    fn tensor_neg() {
        let t1: TensorQuantity<f64, 2, 2, dim::Dimensionless> =
            TensorQuantity::new(nalgebra::SMatrix::from_element(5.0));
        let neg = -t1;
        assert_relative_eq!(neg.matrix[(0, 0)], -5.0);
    }

    #[test]
    fn tensor_mul_raw() {
        let t1: TensorQuantity<f64, 2, 2, dim::Dimensionless> =
            TensorQuantity::new(nalgebra::SMatrix::from_element(3.0));
        let c = t1 * 2.0;
        assert_relative_eq!(c.matrix[(0, 0)], 6.0);
    }

    #[test]
    fn tensor_mul_tensor() {
        let t1: TensorQuantity<f64, 2, 2, dim::Dimensionless> = TensorQuantity::identity();
        let t2: TensorQuantity<f64, 2, 2, dim::Dimensionless> =
            TensorQuantity::new(nalgebra::SMatrix::from_element(7.0));
        let c = t1 * t2;
        assert_relative_eq!(c.matrix[(0, 0)], 7.0);
    }

    #[test]
    fn tensor_default_zero() {
        let t: TensorQuantity<f64, 2, 2, dim::Dimensionless> = TensorQuantity::default();
        assert_relative_eq!(t.matrix[(0, 0)], 0.0);
    }

    #[test]
    fn tensor_into_matrix() {
        let t: TensorQuantity<f64, 2, 2, dim::Dimensionless> = TensorQuantity::identity();
        let raw = t.into_matrix();
        assert_relative_eq!(raw[(0, 0)], 1.0);
    }

    #[test]
    fn tensor_display() {
        let t: TensorQuantity<f64, 2, 2, dim::Force> = TensorQuantity::identity();
        let str_t = format!("{t}");
        assert!(str_t.contains("N"));
    }

    // --- Complex tests ---

    #[test]
    fn complex_quantity() {
        let z = crate::ComplexMeters::new(num_complex::Complex::new(3.0, 4.0));
        assert_relative_eq!(z.value.re, 3.0);
        assert_relative_eq!(z.value.im, 4.0);
    }

    #[test]
    fn complex_velocity() {
        let cv = crate::ComplexFrequency::new(num_complex::Complex::new(1.0, 2.0));
        assert_relative_eq!(cv.value.re, 1.0);
        assert_relative_eq!(cv.value.im, 2.0);
    }

    // --- Prefixed tests ---

    #[test]
    fn prefixed_quantity_types() {
        let _km = Kilometers::new(1.0);
        let _mm = Millimeters::new(1.0);
        let _ms = Milliseconds::new(1.0);
        let _nt = Nanoteslas::new(50.0);
        let _kpa = crate::Kilopascals::new(101.3);
    }

    // --- Dimensionless auto-cast tests ---

    #[test]
    fn dimensionless_into_f64() {
        let ratio: Dimensionless = Dimensionless::new(0.5);
        let val: f64 = ratio.into();
        assert_relative_eq!(val, 0.5);
    }

    #[test]
    fn dimensionless_ref_into_f64() {
        let ratio: Dimensionless = Dimensionless::new(2.71);
        let val: f64 = (&ratio).into();
        assert_relative_eq!(val, 2.71);
    }

    #[test]
    fn dimensionless_from_arithmetic() {
        // m / m = dimensionless
        let ratio: Dimensionless = Meters::new(10.0) / Meters::new(2.0);
        let val: f64 = ratio.into();
        assert_relative_eq!(val, 5.0);
    }

    // --- pow! macro tests ---

    #[test]
    fn si_velocity_from_m_div_s() {
        let vel: Velocity = 5.0 * (m / s);
        assert_relative_eq!(vel.value, 5.0);
    }

    #[test]
    fn si_acceleration() {
        let acc: Acceleration = 9.8 * (m / pow!(s, 2));
        assert_relative_eq!(acc.value, 9.8);
    }

    #[test]
    fn si_area() {
        let area: Area = 3.0 * pow!(m, 2);
        assert_relative_eq!(area.value, 3.0);
    }

    #[test]
    fn si_frequency_from_inverse_s() {
        let freq: Frequency = 1.0 / s;
        assert_relative_eq!(freq.value, 1.0);
    }

    #[test]
    fn si_force() {
        let f: Force = 10.0 * (kg * m / pow!(s, 2));
        assert_relative_eq!(f.value, 10.0);
    }

    #[test]
    fn si_pressure() {
        let p: Pressure = 101_325.0 * (kg / (m * pow!(s, 2)));
        assert_relative_eq!(p.value, 101_325.0);
    }

    #[test]
    fn si_energy() {
        let e: Energy = 4.2 * (kg * pow!(m, 2) / pow!(s, 2));
        assert_relative_eq!(e.value, 4.2);
    }

    #[test]
    fn si_power() {
        let p: Power = 100.0 * (kg * pow!(m, 2) / pow!(s, 3));
        assert_relative_eq!(p.value, 100.0);
    }

    #[test]
    fn si_volume() {
        let vol: Volume = 2.0 * pow!(m, 3);
        assert_relative_eq!(vol.value, 2.0);
    }

    #[test]
    fn si_inverse_meter_wavenumber() {
        let wn: Wavenumber = 1.0 / m;
        assert_relative_eq!(wn.value, 1.0);
    }

    #[test]
    fn si_negative_power() {
        let freq: Frequency = 1.0 * pow!(s, -1);
        assert_relative_eq!(freq.value, 1.0);
    }
}
