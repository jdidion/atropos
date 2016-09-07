from collections import defaultdict
from .report import *
from .seqio import FormatError

def run_serial(reader, modifiers, filters, formatters, writers):
    try:
        n = 0
        total_bp1 = 0
        total_bp2 = 0
        
        for batch_size, batch in reader:
            n += batch_size
            result = defaultdict(lambda: [])
            for record in batch:
                reads, bp = modifiers.modify(record)
                total_bp1 += bp[0]
                total_bp2 += bp[1]
                dest = filters.filter(*reads)
                formatters.format(result, dest, *reads)
            result = dict((path, "".join(strings))
                for path, strings in result.items())
            writers.write_result(result)
        
        return (0, Summary(
            collect_process_statistics(n, total_bp1, total_bp2, modifiers, filters, formatters),
            summarize_adapters(modifiers.get_adapters()),
            modifiers.get_trimmer_classes()
        ).finish())
    
    except KeyboardInterrupt as e:
        logging.getLogger().error("Interrupted")
        return (130, None)
    except IOError as e:
        if e.errno == errno.EPIPE:
            return (1, None)
        raise
    except (FormatError, EOFError) as e:
        logging.getLogger().error("Atropos error", exc_info=True)
        return(1, None)
    finally:
        reader.close()
        writers.close()
