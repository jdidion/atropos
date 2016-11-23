from collections import defaultdict
from atropos.report import *
from atropos.seqio import FormatError

def run_serial(reader, modifiers, filters, formatters, writers):
    report = None
    details = {}
    
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
        
        report = Summary(
            collect_process_statistics(n, total_bp1, total_bp2, modifiers, filters, formatters),
            summarize_adapters(modifiers.get_adapters()),
            modifiers.get_trimmer_classes()
        ).finish()
        
        rc = 0
        details = dict(
            mode='serial',
            threads=1
        )
    
    except KeyboardInterrupt as e:
        logging.getLogger().error("Interrupted")
        rc = 130
    except IOError as e:
        if e.errno == errno.EPIPE:
            rc = 1
        else:
            raise
    except (FormatError, EOFError) as e:
        logging.getLogger().error("Atropos error", exc_info=True)
        rc = 1
    finally:
        reader.close()
        writers.close()
    
    return (rc, report, details)
