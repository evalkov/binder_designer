import os
import psutil
import subprocess

def get_cpu_info():
    # Get the number of physical cores and logical processors (threads)
    cpu_count = os.cpu_count()
    physical_cores = psutil.cpu_count(logical=False)
    return physical_cores, cpu_count

def get_memory_info():
    # Get total and available memory in GB
    mem = psutil.virtual_memory()
    total_memory_gb = mem.total / (1024 ** 3)
    available_memory_gb = mem.available / (1024 ** 3)
    return total_memory_gb, available_memory_gb

def get_disk_io_speed(test_file_path='/tmp/test_io_speed'):
    # Perform a simple disk I/O benchmark using dd
    try:
        write_cmd = f"dd if=/dev/zero of={test_file_path} bs=1G count=1 oflag=dsync"
        read_cmd = f"dd if={test_file_path} of=/dev/null bs=1G"

        print(f"Running write speed test: {write_cmd}")
        write_speed_output = subprocess.check_output(write_cmd, shell=True, stderr=subprocess.STDOUT).decode()
        print(f"Write speed output: {write_speed_output}")
        write_speed = parse_dd_output(write_speed_output)

        print(f"Running read speed test: {read_cmd}")
        read_speed_output = subprocess.check_output(read_cmd, shell=True, stderr=subprocess.STDOUT).decode()
        print(f"Read speed output: {read_speed_output}")
        read_speed = parse_dd_output(read_speed_output)

        return read_speed, write_speed
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e.cmd}")
        print(f"Output: {e.output.decode()}")
        raise
    finally:
        if os.path.exists(test_file_path):
            os.remove(test_file_path)

def parse_dd_output(output):
    # Parse the output of the dd command to find the transfer speed
    for line in output.split('\n'):
        if 'bytes' in line and 'copied' in line:
            parts = line.split(',')
            for part in parts:
                if 'MB/s' in part or 'GB/s' in part:
                    speed_str = part.split()[-2]
                    unit = part.split()[-1]
                    speed = float(speed_str)
                    if unit == 'GB/s':
                        speed *= 1024  # Convert GB/s to MB/s
                    return speed
    raise ValueError("Could not parse dd output for speed")

def calculate_optimal_settings(physical_cores, total_memory_gb, read_speed, write_speed):
    # Set max_workers to the number of physical cores
    max_workers = physical_cores

    # Set batch_size based on available memory, assuming each file is 100MB on average
    file_size_mb = 100
    batch_size = int((total_memory_gb * 1024) * 0.8 / file_size_mb)

    return max_workers, batch_size

def write_to_file(file_path, physical_cores, logical_processors, total_memory_gb, available_memory_gb, read_speed, write_speed, max_workers, batch_size):
    with open(file_path, 'w') as f:
        f.write(f"Physical cores: {physical_cores}\n")
        f.write(f"Logical processors: {logical_processors}\n")
        f.write(f"Total memory: {total_memory_gb:.2f} GB\n")
        f.write(f"Available memory: {available_memory_gb:.2f} GB\n")
        f.write(f"Disk read speed: {read_speed:.2f} MB/s\n")
        f.write(f"Disk write speed: {write_speed:.2f} MB/s\n")
        f.write("\nRecommended settings:\n")
        f.write(f"max_workers: {max_workers}\n")
        f.write(f"batch_size: {batch_size}\n")

def main():
    print("Gathering system information...")

    physical_cores, logical_processors = get_cpu_info()
    total_memory_gb, available_memory_gb = get_memory_info()
    read_speed, write_speed = get_disk_io_speed()

    print(f"Physical cores: {physical_cores}")
    print(f"Logical processors: {logical_processors}")
    print(f"Total memory: {total_memory_gb:.2f} GB")
    print(f"Available memory: {available_memory_gb:.2f} GB")
    print(f"Disk read speed: {read_speed:.2f} MB/s")
    print(f"Disk write speed: {write_speed:.2f} MB/s")

    max_workers, batch_size = calculate_optimal_settings(physical_cores, total_memory_gb, read_speed, write_speed)

    print("\nRecommended settings:")
    print(f"max_workers: {max_workers}")
    print(f"batch_size: {batch_size}")

    write_to_file("compute_settings.txt", physical_cores, logical_processors, total_memory_gb, available_memory_gb, read_speed, write_speed, max_workers, batch_size)

if __name__ == "__main__":
    main()
