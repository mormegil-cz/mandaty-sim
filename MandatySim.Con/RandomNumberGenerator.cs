using System;

namespace MandatySim.Con
{
    public class RandomNumberGenerator : Random
    {
        private const double Multiplier7 = 1.0 / 0x0100_0000_0000_0000L;
        private readonly System.Security.Cryptography.RandomNumberGenerator rng = System.Security.Cryptography.RandomNumberGenerator.Create();

        // private readonly Random random = new();

        /*
        public override int Next(int limit)
        {
            // return random.Next(limit);

            Span<byte> buffer = stackalloc byte[4];
            var max = Int32.MaxValue - ((Int32.MaxValue - (limit - 1)) % limit);
            while (true)
            {
                rng.GetBytes(buffer);
                buffer[3] &= 0x7F;
                var value = BitConverter.ToInt32(buffer);
                if (value <= max) return value % limit;
            }
        }
        */

        protected override double Sample()
        {
            Span<byte> buffer = stackalloc byte[8];
            rng.GetBytes(buffer[..7]);
            return BitConverter.ToInt64(buffer) * Multiplier7;
        }

        public override void NextBytes(byte[] buffer) => rng.GetBytes(buffer);

        public override void NextBytes(Span<byte> buffer) => rng.GetBytes(buffer);
    }
}